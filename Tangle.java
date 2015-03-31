import java.util.*;
import java.io.*;
import java.util.regex.*;

public class Tangle {
  // Pattern for tags
  // Important when replacing them in code.
  // Slashes and dots allowed, but make sure we do not end up
  // doing this within or accross expressions
  // Forbidden: Uppercase Letters, newlines, =, <, >, &, :, ?, !, |, " and '
  public static String half = "[^\\n\\r<>=&?:!|\"'A-Z]+";
  
  // must either be a .java filename or have at least two words.
  public static String TAGPAT=half+" "+half;
  public static String eTAGPAT="<("+TAGPAT+")>";
  
  public static void write(String s, String filename) throws Throwable {
    if (debug)System.out.println("!!! write file "+filename);
    File file = new File(filename);
    if (file.getParentFile() != null) file.getParentFile().mkdirs();
    FileWriter writer = new FileWriter(file);
    
    PrintWriter out = new PrintWriter(filename);
    out.println(s);
    out.close();
  }
  
  public static String[] match(String w, String pat) {
    // null if no match
    Matcher m=Pattern.compile(pat).matcher(w);
    if (!m.find()) return null;
    List<String> allMatches = new ArrayList<String>();
    for (int i = 1; i <= m.groupCount(); i++)
      allMatches.add(m.group(i));
    return allMatches.toArray(new String[0]);
  }
  
  static class Block {
    String name;
    Block(String n, boolean i) {
			ensureBlockWithNameExistsNot(n);
			name = n;
			isFile = i;
			blocks.add(this);
		}
    String content = "";
    boolean isFile = false, used = false;
    
    // expand all blocks found within this block.
    // return false if no blocks are left after this
    boolean expand() {
      if (debug)    System.out.println("!!! --- expand "+name);
      if (debug) System.out.println("!!! before:\n "+content);
      for (Block b : blocks) {
        if (match(content, "<"+b.name+">") != null) {
            b.used = true;
        
            content = content.replaceAll("<"+b.name+">", 
					Matcher.quoteReplacement(b.content));
				// 6.3.2015, 19:56: quoteReplacement required to keep \ (backslash) es
        }
      }
      if (debug) System.out.println("!!! after:\n "+content);
      boolean l = match(content, eTAGPAT) != null;
			
	  if (l) {
	    String[] m = match(content, eTAGPAT);
		blockWithName(m[0]);
		if (debug)System.out.println("!!! <<ERROR>>: "+content.replaceAll(eTAGPAT, "<<ERROR>>"));
	  }
      if (debug) System.out.println("!!! blocks left: "+l);
      
      return l;
    }
  }  
  static List<String> errors = new ArrayList<>();
  static void ensureBlockWithNameExistsNot(String n) {
    for (Block b : blocks) 
      if (b.name.equals(n)) {
        errors.add(line()+"block with name '"+n+"' exists , cannot restart");
        }
  }
  
	static Block blockWithNameOrNull(String n) {
    for (Block b : blocks) 
      if (b.name.equals(n)) return b;
		return null;
	}
	
  static List<String> notFoundErrors = new ArrayList<String>();
  static Block blockWithName(String n) {
        Block b = blockWithNameOrNull(n);
		if (b != null) return b;
        notFoundErrors.add(n);
        return new Block(n, false);
  }
	static Block blockWithNameOrNew(String n, boolean file) {
    Block b = blockWithNameOrNull(n);
		return ( b == null ) ? new Block(n, false) : b;
  }
	
  static List<Block> blocks = new ArrayList<Block>();
  static int linenum = 0;
  static boolean debug, marklines = false;
  static String line() { return "\nLine "+linenum+": "; }
  
  // 0: file to tangle
  // specify any additional argument to enable debugging.
  public static void main(String args[]) throws Throwable {
    BufferedReader br = new BufferedReader(new FileReader(args[0]));
    String line;
    String filename = new File(args[0]).getName();
    Block currentBlock = null;
    
    if (args.length > 1 && args[1].equals("marklines")) {
        marklines = true;
        System.out.println("marklines = true");
    }
    
    debug = (args.length > 1 && !marklines) ||  args.length > 2;
    
    
    List<String> os = new ArrayList<String>();
    while ((line = br.readLine()) != null) {
      linenum++;
      // Process line // only show if debug is on
      if (debug) System.out.println(line);
      
      // Decide what to do
      String[] m;
      
	  // error (typo): normal line starting with [ (also filter out < if only a single word)
      if ((m = match(line, "^\\[([^\\]]*)\\]=[ \t]*$")) != null ||
		(m = match(line, "^<("+TAGPAT+")>(\\+?)=[ \t]*$")) != null) {
        System.err.println(line()+"'"+line+"'\n\terror: normal line starting with tag, must be indented");
		System.exit(1);
      }
	  
      // <>= // or +=
      if ((m = match(line, "^[ \t]+<("+TAGPAT+")>(\\+?)=[ \t]*$")) != null) {
        // Start <>= block
        if (m[1].equals("+")) {
          if (debug)System.out.println("!!! Discovered <"+m[0]+">+= line");
          currentBlock = blockWithNameOrNew(m[0], false); // 6.3.: Create blocks with +=.
					// Even when not previously used with =
				}
        else {
          if (debug)System.out.println("!!! Discovered <"+m[0]+">= line");
          currentBlock = new Block(m[0], false);
        }
        continue;
      }
      else if ((m = match(line, "^[ \t]+<("+TAGPAT+")>\\+")) != null) {
        System.out.println(line() +", '"+line+"', looks like a tag definition <> but is missing = from +=. Add a space before + if you did not mean that.");
        System.exit(1);
      }
      
      // []=
      if ((m = match(line, "^[ \t]+\\[([^\\]]*)\\]=[ \t]*$")) != null) {
        // Start <>= block
        if (debug)System.out.println("!!! Discovered ["+m[0]+"]= line");
        currentBlock = new Block(m[0], true);
        continue;
      }
      else  if ((m = match(line, "^[ \t]+\\[([^\\]]*)\\]")) != null) {
        System.out.println(line() +", '"+line+"', begins like a file definition [] but is missing =");
        System.exit(1);
      }
      
      // $
      if ((m = match(line, "^[ \t]+\\$(.*)$")) != null) {
        System.out.println("!!! OS $: "+m[0]);
        //java.lang.Runtime.getRuntime().exec(m[0]);
        os.add(m[0]);
        continue;
      }
      
      // End of block
      if (match(line, "^[^ \t]") != null && currentBlock != null) {
        if (debug)System.out.println("!!! end of block "+currentBlock.name);
        currentBlock = null;
        continue;
      }
      
      // Add to block
      if (currentBlock != null) {
        if (marklines) 
            currentBlock.content += "#line "+linenum+" \""+filename+"\"\r\n";
       currentBlock.content += line+"\r\n";
      }
      // eol
    }
    // eof
    br.close();
    
    // Debug
    System.out.println("!!! === Blocks === ");
		if (debug)
    for (Block b : blocks) {
      System.out.println("!!! --- Block '"+b.name+"' --- ");
      System.out.println(b.content);
    }
		else
			System.out.println("<omitted>");
    
    System.out.println("!!! === end of Blocks === ");
    
    // Expand
    List<Block> pb = new ArrayList<Block>();
    pb.addAll(blocks);
    while (!pb.isEmpty()) {
      List<Block> pbn = new ArrayList<Block>();
      for (Block b : pb) {
        if (b.expand()) pbn.add(b);
      }
      pb = pbn;
    }

    // Nothing left to expand. 
    
    // Errors?
    
    // Undefined blocks (failed replacements)
    boolean ok = true;
    if (!notFoundErrors.isEmpty()) {
        for (String s : notFoundErrors) 
            System.err.println("block with name '"+s+"' not found among existing, could not do all replacements");
        ok = false;
    }
    
    
    if (!errors.isEmpty()) {
        for (String s : errors) 
            System.err.println(s);
        ok = false;
    }
    
    // Unused non-file blocks
    for (Block b : blocks) 
        if (!b.isFile && !b.used) {
            System.err.println("block with name '"+b.name+"' never used");
            ok = false;
        }
    if (!ok) System.exit(1);
        
    // Write
    for (Block b : blocks) 
      if (b.isFile)
        write(b.content, b.name);
    
    // os
    for (String o : os) {
      System.out.println("!!! OS run $: "+o);
      //java.lang.Runtime.getRuntime().exec(o).waitFor(); // must wait // cannot get output with this?
			ProcessBuilder pbx = new ProcessBuilder(o.split(" ")).inheritIO().redirectErrorStream(true);
			pbx.redirectOutput(ProcessBuilder.Redirect.INHERIT);
			pbx.start().waitFor();
    }
    
  }
}
