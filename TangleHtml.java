import java.util.*;
import java.io.*;
import java.util.regex.*;
import java.nio.charset.*;
import java.nio.file.*;

public class TangleHtml {
  // Pattern for tags
  // Important when replacing them in code.
  // Slashes and dots allowed, but make sure we do not end up
  // doing this within or accross expressions
  public static String TAGPAT="[^<>=&|\"']+";
  public static String eTAGPAT="<"+TAGPAT+">";
  
  public static void write(String s, String filename) throws Throwable {
    System.out.println("!!! write file "+filename);
    File file = new File(filename);
    if (file.getParentFile() != null) file.getParentFile().mkdirs();
    //FileWriter writer = new FileWriter(file);
    
    Writer out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "UTF-8"));
	try {
		out.write(s);
	} finally {
		out.close();
	}
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
  private static String OS = null;
  public static String getOsName()
  {
      if(OS == null) { OS = System.getProperty("os.name"); }
      return OS;
  }
  public static boolean isWindows()
  {
      return getOsName().startsWith("Windows");
  }
  static String readResourceFile(String path) throws IOException {
		java.net.URL main = TangleHtml.class.getResource(path);
		String fn = main.getPath();

        // on windows
        if (isWindows()) fn = fn.substring(1);
		
		//byte[] encoded = Files.readAllBytes(Paths.get(main.getURI()));
		
		System.out.println("|"+fn);
		return //new Scanner(new File(
			new String(Files.readAllBytes(Paths.get(fn)), "UTF-8");
		//)).useDelimiter("\\Z").next();
		//new String(encoded, Charset.defaultCharset());
	}
  
  static int linenum = 0;
  
  public static void main(String args[]) throws Throwable {
    BufferedReader br = new BufferedReader(
           new InputStreamReader(
                      new FileInputStream(args[0]), "UTF8"));;
    String line;
    String out = readResourceFile("TangleHtml.html");
    boolean b = false; // in block
    while ((line = br.readLine()) != null) {
      linenum++;
      // Process line
      
      if (
				// <>= // or +=
				match(line, "^[ \t]+<("+TAGPAT+")>(\\+?)=[ \t]*$") != null
				// []=
				|| match(line, "^[ \t]+\\[("+TAGPAT+")\\]=[ \t]*$") != null
				// $
				|| match(line, "^[ \t]+\\$") != null
				) {
				if (b)  out += "</pre>\r\n";
        out += "<pre>\r\n";b = true;
      }
      
			if (!b &&
				// or just any indented block
				match(line, "^[ \t]+") != null
				) {
        out += "<pre>\r\n";b = true;
      }
			
      // End of block
      if (match(line, "^[^ \t]") != null && b) {
        out += "</pre>\r\n"; b = false;
      }
			
			if (b) out += line.replaceAll("<", "&lt;").replaceAll(">", "&gt;");
			else out += line;
			out += "\r\n";
		}
    // eof
    br.close();
    
    write(out, args[0].replaceAll("\\..*", ".html"));
  }
}
