<h2>Sampler</h2>
Samplers make random samples, which are used for Monte Carlo rendering. The 
samples always need to lie in the range [0,1]. Various versions such 
as purely random, jittered, or low discrepancy samplers could be 
implemented.
<p>
makeSamples makes an array of samples. The samples need to lie in the range [0,1]^d,
where d is the dimensionality of the samples.
The number of returned samples may differ from the number of desired samples n.
	<make samples>=
	public float[][] makeSamples(int n, int d);  
	[rt/Sampler.java]= 
	package rt;
	public interface Sampler {
		<make samples>
	}
		
	[rt/SamplerFactory.java]= 
	package rt;

	public interface SamplerFactory {
		public Sampler make();
	}
	
<h3>OneSampler</h3>
Returns always one sample at 0.5 in all dimensions.
	[rt/samplers/OneSampler.java]= 
	package rt.samplers;

	<common imports>

	public class OneSampler implements Sampler {		
		public float[][] makeSamples(int n, int d)
		{
			float[][] samples = new float[1][d];
			for(int i=0; i<d; i++)
				samples[0][i] = 0.5f;
			
			return samples;
		}
	}
		
	[rt/samplers/OneSamplerFactory.java]= 
	package rt.samplers;

	<common imports>
	public class OneSamplerFactory implements SamplerFactory {
		public Sampler make() {return new OneSampler();}
	}
		
<h3>RandomSampler</h3>
Makes uniform random samples in the range [0,1].
	[rt/samplers/RandomSampler.java]= 
	package rt.samplers;

	<common imports>
	public class RandomSampler implements Sampler {

		Random random = new Random();
		public float[][] makeSamples(int n, int d)
		{
			float samples[][] = new float[n][d];
			
			for(int i=0; i<n; i++) for(int j=0; j<d; j++)
			{
				samples[i][j] = random.nextFloat();
			}
			return samples;
		}
		
	}
		
	[rt/samplers/RandomSamplerFactory.java]= 
	package rt.samplers;
	<common imports>
	public class RandomSamplerFactory implements SamplerFactory {
		public Sampler make() {return new RandomSampler();}
	}
	
    
<h2>Film</h2>
This and the next section are concerned with the pseudocode
    set pixel color
of the basic raytracing algorithm.
<p>
A film captures the samples generated for (sub) pixels
	<add sample>=
	public void addSample(double x, double y, Spectrum s);
and accretes them in discrete pixels, much like a photo sensor in a digital camera integrates
the light arriving at the square shaped sensors.
The final image is a rectangular array of light spectra.
	<retrieve array>=
	public Spectrum[][] getImage();
	public int getWidth();
	public int getHeight();
	
So a film stores a 2D grid of Spectrum representing an image.
Rendered samples can be added one by one to a film. Samples are
filtered using some filter (depending on the implementation of this 
interface) when added.
	[rt/Film.java]= 
	package rt;
	public interface Film {
		<add sample>
		<retrieve array>
	}
	
<h3>BoxFilterFilm</h3>	
This film uses a box filter when accumulating samples on a film.
A box filter means that samples s contribute only to the pixel (x,y) that they lie in.
	<add sample to film>=
	int idx_x = (int)x;
	int idx_y = (int)y;
    
	unnormalized[idx_x][idx_y].add(s);
	nSamples[idx_x][idx_y]++;
    
Sample values are simply averaged for the final image.
	<average sample values>=
	image[i][j] = new Spectrum(unnormalized[i][j]).mult(1.f/nSamples[i][j]);

	[rt/films/BoxFilterFilm.java]= 
	package rt.films;
    <common imports>
	public class BoxFilterFilm implements Film {
		private Spectrum[][] image;
		public int width, height;
		private Spectrum[][] unnormalized;
		private float nSamples[][];
		
		public BoxFilterFilm(int width, int height)
		{
			this.width = width;
			this.height = height;
			image = new Spectrum[width][height];
			unnormalized = new Spectrum[width][height];
			nSamples = new float[width][height];
			
			for(int i=0; i<width; i++) for(int j=0; j<height; j++) {
                image[i][j] = new Spectrum();
                unnormalized[i][j] = new Spectrum();
            }
		}
		
		public void addSample(double x, double y, Spectrum s)
		{
			if((int)x>=0 && (int)x<width && (int)y>=0 && (int)y<height) {
				<add sample to film>
			}
		}
		
		public int getWidth() {return width;}
		public int getHeight() {return height;}
		
		public Spectrum[][] getImage()
		{
			for(int i=0; i<width; i++) for(int j=0; j<height; j++) {
                <average sample values>
            }
			return image;
		}
	}
		
<h2>Tonemapper</h2>
A Tonemapper compresses a raw rendered Film to an image that can be displayed on typical 8-bit displays.
	[rt/Tonemapper.java]= 
	package rt;
	<common imports>

	public interface Tonemapper {
		BufferedImage process(Film film);
	}
	
Tone maps a film by clamping all color channels to range [0,1].
	<do clamping>=
	Spectrum s = filmImg[i][j].clamp();
    
Then we convert these to 24 bit (24 bpp) 8-bit per channel colors.
    <convert spectrum to rgb integer>=
     ((int)(255.f*s.r) << 16) 
    | ((int)(255.f*s.g) << 8) 
    | ((int)(255.f*s.b))
    
	[rt/tonemappers/ClampTonemapper.java]= 
	package rt.tonemappers;
	
	<common imports>

	public class ClampTonemapper implements Tonemapper {
		public BufferedImage process(Film film)
		{
			BufferedImage img = new BufferedImage(film.getWidth(), film.getHeight(), BufferedImage.TYPE_3BYTE_BGR);
			
			Spectrum[][] filmImg = film.getImage();
			for(int i=0; i<film.getWidth(); i++) for(int j=0; j<film.getHeight(); j++) {
				<do clamping>
                
				img.setRGB(i, film.getHeight()-1-j, 
                    <convert spectrum to rgb integer>
				);
			}
			return img;
		}
	}

<h2>Test Scenes</h2>
<h3>Rendered Scenes</h3>

The following class encapsulates some defaults shared by most scenes.
	[rt/basicscenes/AbstractScene.java]= 
	package rt.basicscenes;
	<common imports>
	public abstract class AbstractScene extends Scene {
		public void setDimensions(int w, int h) {
			width = w; height = h;
			film = new BoxFilterFilm(width, height);
			camera = new FixedCamera(width, height);
		}
		
		public void setDimensions(int wh) {
			setDimensions(wh, wh);
		}
		
		public void setSPP(int s) {
			SPP = s;
			outputFilename = new String("output/"+this.getClass().getName());
			outputFilename = outputFilename + " " + String.format("%d", SPP) + "SPP";
			if (s > 1) samplerFactory = new RandomSamplerFactory();
			else samplerFactory = new OneSamplerFactory();
		}
		
		public AbstractScene() {
			setDimensions(512);
			setSPP(1);
			tonemapper = new ClampTonemapper();
			samplerFactory = new OneSamplerFactory();
			integratorFactory = new PointLightIntegratorFactory();
			lightList = new LightList();
			LightGeometry pointLight = new PointLight(new Vector3f(0.f, 0.f, 3.f), new Spectrum(15.f));
			lightList.add(pointLight);
		}
	}
		
Here are two tests demonstrating that the FixedCamera works:
	[rt/basicscenes/Box.java]= 
	package rt.basicscenes;
	<common imports>
	public class Box extends AbstractScene {
		public Box() {
			root = new CSGCappedZTunnel();
			<box light>
		}
	}
	<test scenes>+=
	new Box(),
	<box light>=
	LightGeometry pointLight = new PointLight(new Vector3f(0.f, 0.f, 3.f), new Spectrum(10.f, 10.f, 10.f));
	lightList = new LightList();
	lightList.add(pointLight);
		
<img src="output/rt.basicscenes.Dodecahedron 1SPP.png"></img>
	[rt/basicscenes/Dodecahedron.java]= 
	package rt.basicscenes;
	<common imports>
	public class Dodecahedron extends AbstractScene {
		public Dodecahedron()
		{
			root = CSGNode.add(
					new CSGPlane(new Vector3f(0.f, 1.f, 0.f), 1.f),
					new CSGPlane(new Vector3f(0.f, 0.f, 1.f), 1.f),
					new CSGDodecahedron()
				);
		}
	}	
	<test scenes>+=
	new Dodecahedron(),
	
	[rt/testscenes/PinholeCameraScene.java]= 
	package rt.testscenes;
	<common imports>
	public abstract class PinholeCameraScene extends AbstractScene {
		public Vector3f eye, lookAt, up;
		
		public void setDimensions(int w, int h) {
			super.setDimensions(w, h);
			setCamera(eye == null ? new Vector3f() : eye, 
					lookAt == null ? new Vector3f() : lookAt, 
							up == null ? new Vector3f() : up);
		}
		
		public PinholeCameraScene() {}
		public PinholeCameraScene(Vector3f eye) {setCamera(eye, new Vector3f(), <up vector>);}
		
		public void setCamera(Vector3f eye,
				Vector3f lookAt,
				Vector3f up) {
			this.eye = eye;
			this.lookAt = lookAt;
			this.up = up;
			float verticalFovInDegrees = 60.f;
			float aspect = (width * 1.f)/height;
			camera = new PinholeCamera(eye, lookAt, up, verticalFovInDegrees, aspect, width, height);
		}
	}
			
A test baseclass for objects that looks at the origin and uses the debug integrator
	[rt/testscenes/ObjectTest.java]= 
	package rt.testscenes;
	<common imports>
	public abstract class ObjectTest extends PinholeCameraScene {
		public ObjectTest(Vector3f eye) {
			setCamera(eye, new Vector3f(), <up vector>);
			integratorFactory = new DebugIntegratorFactory();
		}
	}
Same thing but with the normal debug integrator.
	[rt/testscenes/ObjectNormalTest.java]= 
	package rt.testscenes;
	<common imports>
	public abstract class ObjectNormalTest extends ObjectTest {
		public ObjectNormalTest(Vector3f eye) {
			super(eye);
			integratorFactory = new NormalDebugIntegratorFactory();
		}
	}
Same thing but with the uv debug integrator.
	[rt/testscenes/ObjectUVTest.java]= 
	package rt.testscenes;
	<common imports>
	public abstract class ObjectUVTest extends ObjectTest {
		public ObjectUVTest(Vector3f eye) {
			super(eye);
			integratorFactory = new UVDebugIntegratorFactory();
		}
	}
Test scene for pinhole camera specifications.
<img src="output/rt.testscenes.CameraTestScene 1SPP.png"></img>
	[rt/testscenes/CameraTestScene.java]= 
	package rt.testscenes;
	<common imports>
	public class CameraTestScene extends PinholeCameraScene {
		public CameraTestScene()
		{
			setDimensions(1280, 720);
			Vector3f eye = new Vector3f(0.5f, 0.5f, 3.f);
			Vector3f lookAt = new Vector3f(0.5f, 0.f, 0.f);
			Vector3f up = new Vector3f(0.2f, 1.f, 0.f);
			setCamera(eye, lookAt, up);
			root = new CSGCappedZTunnel();
			<box light>
		}
	}
	<test scenes>+=
	new CameraTestScene(),
	
<h3>Beautiful Scenes</h3>
<h4>InstancingTeapots</h4>
<img src="output/rt.testscenes.InstancingTeapots 1SPP.png"></img>
	[rt/testscenes/InstancingTeapots.java]= 
	package rt.testscenes;
	<common imports>
	public class InstancingTeapots extends PinholeCameraScene {
		public InstancingTeapots()
		{	
			super(new Vector3f(0.f,0.f,2.f));
			setDimensions(256);
			
			// List of objects
			IntersectableList objects = new IntersectableList();	
					
			// Box
			CSGPlane plane = new CSGPlane(new Vector3f(0.f, 1.f, 0.f), 1.f);
			plane.material = new Diffuse(new Spectrum(0.f, 0.8f, 0.8f));
			objects.add(plane);		
			
			plane = new CSGPlane(new Vector3f(0.f, 0.f, 1.f), 1.f);
			plane.material = new Diffuse(new Spectrum(0.3f, 0.8f, 0.8f));
			objects.add(plane);
			
			plane = new CSGPlane(new Vector3f(-1.f, 0.f, 0.f), 1.f);
			plane.material = new Diffuse(new Spectrum(1.f, 0.8f, 0.8f));
			objects.add(plane);
			
			plane = new CSGPlane(new Vector3f(1.f, 0.f, 0.f), 1.f);
			plane.material = new Diffuse(new Spectrum(0.f, 0.8f, 0.0f));
			objects.add(plane);
			
			plane = new CSGPlane(new Vector3f(0.f, -1.f, 0.f), 1.f);
			plane.material = new Diffuse(new Spectrum(0.8f, 0.8f, 0.8f));
			objects.add(plane);
			
			// Add objects
			Mesh mesh = ObjReader.read("obj/teapot.obj", 1.f);
			Matrix4f t = new Matrix4f();
			t.setIdentity();
			
			// Instance one
			t.setScale(0.5f);
			t.setTranslation(new Vector3f(0.f, -0.35f, 0.f));
			Instance instance = new Instance(mesh, t);
			objects.add(instance);	
			
			// Instance two
			t.setScale(0.5f);
			t.setTranslation(new Vector3f(0.f, 0.25f, 0.f));
			Matrix4f rot = new Matrix4f();
			rot.setIdentity();
			rot.rotX((float)Math.toRadians(30.f));
			t.mul(rot);
			instance = new Instance(mesh, t);
			objects.add(instance);
					
			root = objects;
			
			// List of lights
			lightList = new LightList();
			
			LightGeometry light = new PointLight(new Vector3f(0.f,0.8f,0.8f), new Spectrum(3.f, 3.f, 3.f));
			lightList.add(light);
			
			light = new PointLight(new Vector3f(-0.8f,0.2f,1.f), new Spectrum(1.5f, 1.5f, 1.5f));
			lightList.add(light);		
		}
	}
	<beautiful scenes>+=
	new InstancingTeapots(),
	
<h4>Blinn</h4>
Simple scene using a Blinn material.	
<img src="output/rt.testscenes.BlinnTest 1SPP.png"></img>
	[rt/testscenes/BlinnTest.java]= 
	package rt.testscenes;
	<common imports>
	public class BlinnTest extends PinholeCameraScene {

		public BlinnTest()
		{
			setDimensions(512,512);
			setSPP(1);
			Vector3f eye = new Vector3f(0.f, 0.f, 3.f);
			Vector3f lookAt = new Vector3f(0.f, 0.f, 0.f);
			setCamera(eye, lookAt, <up vector>);
			
			// Specify which integrator and sampler to use
			integratorFactory = new PointLightIntegratorFactory();
			samplerFactory = new OneSamplerFactory();

			// Ground plane
			CSGPlane groundPlane = new CSGPlane(new Vector3f(0.f, 1.f, 0.f), 1.f);
			
			// Sphere with Blinn material
			CSGSphere sphere = new CSGSphere();
			sphere.material = new Blinn(new Spectrum(.8f, 0.f, 0.f), new Spectrum(.4f, .4f, .4f), 50.f);
			
			root = new IntersectableList().add(groundPlane, sphere);
			
			// Light sources
			LightGeometry pl1 = new PointLight(new Vector3f(.5f, .5f, 2.f), new Spectrum(1.f, 1.f, 1.f));
			LightGeometry pl2 = new PointLight(new Vector3f(-.75f, .75f, 2.f), new Spectrum(1.f, 1.f, 1.f));
			lightList = new LightList();
			lightList.add(pl1);
			lightList.add(pl2);
		}
	}
	<beautiful scenes>+=
	new BlinnTest(),
   
    
<h4>AreaLightScene</h4>
<img src="output/rt.testscenes.AreaLightScene 32SPP.png"></img>
	[rt/testscenes/AreaLightScene.java]= 
	package rt.testscenes;
	<common imports>
	public class AreaLightScene extends PinholeCameraScene {
		public AreaLightScene()
		{
			super(new Vector3f(0.f, 0.f, 5.f));
			setDimensions(512);
            setSPP(32);
			integratorFactory = new AreaLightIntegratorFactory();
			
            // Arealight requires random samples
            samplerFactory = new RandomSamplerFactory();
            
            <ground and back plane>
			
			root =  new IntersectableList().add(
				groundPlane,
				backPlane,
				new CSGSphere()
				);
			
			lightList = new LightList();
			lightList.add(new AreaLight(
                new Vector3f(0,3.f,0), 
                new Vector3f(4.f,0,0), 
                new Vector3f(0,4.f,0), 
                new Spectrum(80.f)));
		}
	}
	<beautiful scenes>+=
	new AreaLightScene(),
    
<h4>ShadowScene</h4>
In the following scenes, we have a backplane at z=
	<backplane z>=
	-3.15f
	
and a groundplane at y
	<groundplane y>=
	-1.5f
	
The sphere in the center is at the origin and the camera at 0,0,5.

The light sources are at y = 3, z = 2 and shifted left and right:
	<some light sources>=
	LightGeometry pointLight1 = new PointLight(
		new Vector3f(-1.f, 3.f, 2.f), 
		new Spectrum(44.f));
	LightGeometry pointLight2 = new PointLight(
		new Vector3f(1.f, 3.f, 2.f),
		new Spectrum(44.f));
	lightList = new LightList();
	lightList.add(pointLight1);
	lightList.add(pointLight2);
	
<img src="output/rt.testscenes.ShadowScene 1SPP.png"></img>
	[rt/testscenes/ShadowScene.java]= 
	package rt.testscenes;
	<common imports>
	public class ShadowScene extends PinholeCameraScene {
		public ShadowScene()
		{
			super(new Vector3f(0.f, 0.f, 5.f));
			setDimensions(512);
			integratorFactory = new WhittedIntegratorFactory();
			<ground and back plane>
			
			root =  new IntersectableList().add(
				groundPlane,
				backPlane,
				new CSGSphere()
				);
			
			<some light sources>
		}
	}
	<beautiful scenes>+=
	new ShadowScene(),
	
	<ground and back plane>=
	XYZGrid grid = new XYZGrid(new Spectrum(0.2f, 0.f, 0.f), new Spectrum(1.f, 1.f, 1.f), 0.1f, new Vector3f(0.f, 0.3f, 0.f));
	CSGPlane groundPlane = new CSGPlane(new Vector3f(0.f, 1.f, 0.f), - <groundplane y>);
	groundPlane.material = grid;
	CSGPlane backPlane = new CSGPlane(new Vector3f(0.f, 0.f, 1.f), - <backplane z>);
	backPlane.material = grid;		

<img src="output/rt.testscenes.RefractionScene 1SPP.png"></img>
	[rt/testscenes/RefractionScene.java]= 
	package rt.testscenes;
	<common imports>
	public class RefractionScene extends PinholeCameraScene {
		public RefractionScene()
		{
			super(new Vector3f(0.f, 0.f, 5.f));
			setDimensions(512);
			integratorFactory = new WhittedIntegratorFactory();
			<ground and back plane>
			
			root =  new IntersectableList().add(
				groundPlane,
				backPlane,
				new CSGSphere(new Vector3f(), 1, new RefractiveOnly(1.3f))
				);
			
			<some light sources>
		}
	}
	<beautiful scenes>+=
	new RefractionScene(),
	
<img src="output/rt.testscenes.RefractionScene2 1SPP.png"></img>
	[rt/testscenes/RefractionScene2.java]= 
	package rt.testscenes;
	<common imports>
	public class RefractionScene2 extends PinholeCameraScene {
		public RefractionScene2()
		{
			super(new Vector3f(0.f, 0.f, 5.f));
			setDimensions(512);
			integratorFactory = new WhittedIntegratorFactory();
			<ground and back plane>
			
			root =  new IntersectableList().add(
				groundPlane,
				backPlane,
				CSGNode.intersect(
					new CSGPlane(new Vector3f(0,0,1.f),0,new RefractiveOnly(1.3f)),
					new CSGPlane(new Vector3f(1.f,0,0),0,new RefractiveOnly(1.3f))
					)
				);
			
			<some light sources>
		}
	}
	<beautiful scenes>+=
	new RefractionScene2(),
	
<h4>ReflectionScene</h4>
<img src="output/rt.testscenes.ReflectionScene 1SPP.png"></img>
	[rt/testscenes/ReflectionScene.java]= 
	package rt.testscenes;
	<common imports>
	public class ReflectionScene extends PinholeCameraScene {
		public ReflectionScene() {
			super(new Vector3f(0.f, 0.f, 5.f));
			setDimensions(512);
			integratorFactory = new WhittedIntegratorFactory();
			<ground and back plane>

			root =  new IntersectableList().add(
				groundPlane,
				backPlane,
				new CSGSphere(new Vector3f(), 1, new Reflective()),
				new CSGSphere(new Vector3f(2.f, 0, 0), 1)
				);
			
			<some light sources>
		}
	}
	<beautiful scenes>+=
	new ReflectionScene(),
	
<h4>CSGScene</h4>
<img src="output/rt.testscenes.CSGScene 1SPP.png"></img>
	[rt/testscenes/CSGScene.java]= 
	package rt.testscenes;
	<common imports>
	public class CSGScene extends PinholeCameraScene {
		public CSGScene()
		{
			setSPP(32);
			setDimensions(640, 360);

			Vector3f eye = new Vector3f(0.f, 0.f, 5.f);
			Vector3f lookAt = new Vector3f(0.f, -.5f, 0.f);
			setCamera(eye, lookAt, <up vector>);
			
			integratorFactory = new WhittedIntegratorFactory();
			
			Material refractive = new Refractive(1.3f);
			
			// Make a conical "bowl" by subtracting cross-sections of two cones
			CSGSolid outerCone = coneCrossSection(60.f, refractive);
			// Make an inner cone and subtract it
			Matrix4f trafo = new Matrix4f();
			trafo.setIdentity();
			trafo.setTranslation(new Vector3f(0.f, 0.f, 0.25f));
			CSGInstance innerCone = new CSGInstance(outerCone, trafo);		
			CSGSolid doubleCone = new CSGNode(outerCone, innerCone, CSGNode.OperationType.SUBTRACT);
			
			// Place it in the scene
			Matrix4f rot = new Matrix4f();
			rot.setIdentity();
			rot.rotX(-(float)Math.PI/2.f);
			Matrix4f trans = new Matrix4f();
			trans.setIdentity();
			trans.setTranslation(new Vector3f(-1.5f, -1.5f, 0.f));
			trans.mul(rot);		
			doubleCone = new CSGInstance(doubleCone, trans);
			
			// Something like a"soap bar"
			Material yellow = new Diffuse(new Spectrum(1.f, 0.8f, 0.2f));
			CSGSolid soap = new CSGUnitCylinder(yellow);
			CSGSolid cap = new CSGTwoSidedInfiniteCone(yellow);
			// Smoothen the edges
			trans.setIdentity();
			trans.m23 = -0.8f;
			CSGSolid cap1 = new CSGInstance(cap, trans); 
			soap = new CSGNode(soap, cap1, CSGNode.OperationType.INTERSECT);
			trans.m23 = 1.8f;
			CSGSolid cap2 = new CSGInstance(cap, trans); 
			soap = new CSGNode(soap, cap2, CSGNode.OperationType.INTERSECT);
			
			// Transform it and place it in the scene
			Matrix4f scale = new Matrix4f();
			// Make it elliptical and rotate a bit around the cylinder axis
			scale.setIdentity();
			scale.m11 = 0.5f;
			scale.m22 = 0.5f;
			trafo = new Matrix4f();
			trafo.rotZ((float)Math.toRadians(-20));
			trafo.mul(scale);
			// Rotate it "up"
			rot = new Matrix4f();
			rot.setIdentity();
			rot.rotX(-(float)Math.PI/2.f);		
			rot.mul(trafo);
			// Place in scene by translating
			trans = new Matrix4f();
			trans.setIdentity();
			trans.setTranslation(new Vector3f(1.5f, -1.5f, 1.f));
			trans.mul(rot);
			soap = new CSGInstance(soap, trans);
			
			// Ground and back plane
			XYZGrid grid = new XYZGrid(new Spectrum(0.2f, 0.f, 0.f), new Spectrum(1.f, 1.f, 1.f), 0.1f, new Vector3f(0.f, 0.3f, 0.f));
			CSGPlane groundPlane = new CSGPlane(new Vector3f(0.f, 1.f, 0.f), 1.5f);
			groundPlane.material = grid;
			CSGPlane backPlane = new CSGPlane(new Vector3f(0.f, 0.f, 1.f), 3.15f);
			backPlane.material = grid;		
			
			// Sphere
			CSGSphere sph = new CSGSphere();
			sph.material = new Refractive(1.01f);
			
			// Set the root node for the scene
			root =  new IntersectableList().add(
				doubleCone,
				soap,
				groundPlane,
				backPlane,
				sph
				);
			
			<some light sources>
		}
		
		<cone cross section>
	}
	<beautiful scenes>+=
	new CSGScene(),
	
Where
	<cone cross section>=
	private CSGSolid coneCrossSection(float a, Material material)
	{
		// Makes a two-sided infinite cone with apex angle 90 degrees
		CSGTwoSidedInfiniteCone doubleCone = new CSGTwoSidedInfiniteCone(material);
		// Scaling factor along the cone axis corresponding to apex angle
		float s = (float)Math.tan((90-a/2)/180.f*(float)Math.PI);
		
		// Scale and translate cone
		Matrix4f scale = new Matrix4f();
		scale.setIdentity();
		scale.m22 = s;
		Matrix4f trans = new Matrix4f();
		trans.setIdentity();
		trans.setTranslation(new Vector3f(0.f, 0.f, -s));
		trans.mul(scale);
		CSGInstance scaledCone = new CSGInstance(doubleCone, trans);
		
		// Cut off at z=0 and z=1
		CSGNode out = new CSGNode(scaledCone, new CSGPlane(new Vector3f(0.f, 0.f, -1.f), 0.f, material), CSGNode.OperationType.INTERSECT);
		out = new CSGNode(out, new CSGPlane(new Vector3f(0.f, 0.f, 1.f), -1.f, material), CSGNode.OperationType.INTERSECT);
		
		return out;
	}
Makes a "horizontal" cross section through a cone with apex angle a (in degrees).
The bottom plane is at z=0, the top at z=1. The radius of the bottom circle 
in the cross section is one (the top circle is bigger depending on the apex angle).
<img src=ccs.jpg></img>

<h4>Refractive Glass Block</h4>
Test scene for refractive glass block.
Features a background plane with a debug world position material:
	<construct backplane>=
	CSGPlane backPlane = new CSGPlane(new Vector3f(0.f, 1.f, 0.f), 0.f);
	backPlane.material = new DebugMaterial();
and a glass block occupying half of it (notice that plane-origin offsets are given in the direction opposite to the normal)
	<construct glassblock>=
	Material refractive = new Refractive(1.3f);
	CSGSolid glassblock = CSGNode.intersect(
		new CSGPlane(new Vector3f(0.f, 1.f, 0.f), -1.f, refractive),
		new CSGPlane(new Vector3f(0.f, -1.f, 0.f), 0.5f, refractive),
		new CSGPlane(new Vector3f(1.f, 0.f, 0.f), 0.f, refractive)
		);
	
<img src=glassplane.JPG></img>
<img src="output/rt.testscenes.RefractiveGlassplane 1SPP.png"></img>	
	[rt/testscenes/RefractiveGlassplane.java]= 
	package rt.testscenes;
	<common imports>
	public class RefractiveGlassplane extends PinholeCameraScene {
		public RefractiveGlassplane()
		{
			super(new Vector3f(0.f, 3.f, 3.f));
			setDimensions(512);
			
			integratorFactory = new WhittedIntegratorFactory();
			
			<construct backplane>
			<construct glassblock>
			
			root = new IntersectableList().add(
				backPlane,
				glassblock
				);
			
			lightList = new LightList();
			lightList.add(
				new PointLight(new Vector3f(3.f, 3.f, 0.f), new Spectrum(14.f))
			);
		}
		
	}
	<beautiful scenes>+=
	new RefractiveGlassplane(),
	
<h4>Refractive Sphere</h4>
Test scene for refractive objects, renders a sphere in front of a planar background.
<img src="output/rt.testscenes.RefractiveSphere 32SPP.png"></img>	
	[rt/testscenes/RefractiveSphere.java]= 
	package rt.testscenes;
	<common imports>
	public class RefractiveSphere extends PinholeCameraScene {
			
		public RefractiveSphere()
		{
			setDimensions(512);
			setSPP(32);
			
			// Specify which camera, film, and tonemapper to use
			Vector3f eye = new Vector3f(0.f, 0.f, 3.f);
			Vector3f lookAt = new Vector3f(0.f, 0.f, 0.f);
			setCamera(eye, lookAt, <up vector>);
			
			// Specify which integrator and sampler to use
			integratorFactory = new WhittedIntegratorFactory();
			samplerFactory = new RandomSamplerFactory();		
			
			Material refractive = new Refractive(1.3f);

			
			// Ground and back plane
			// A grid with red and white lines, 
			// line thickness 0.01, zero offset shift, and tile size 0.125, all in world coordinates
			//XYZGrid grid = new XYZGrid(new Spectrum(0.2f, 0.f, 0.f), new Spectrum(1.f, 1.f, 1.f), 0.01f, new Vector3f(0.f, 0.f, 0.f), 0.125f);
			XYZGrid grid = new XYZGrid(new Spectrum(0.2f, 0.f, 0.f), new Spectrum(1.f, 1.f, 1.f), 0.1f, new Vector3f(0.f, 0.3f, 0.f));
			//CSGPlane backPlane = new CSGPlane(new Vector3f(0.f, 0.f, 1.f), 2.15f);
			CSGPlane backPlane = new CSGPlane(new Vector3f(0.f, 0.f, 1.f), 3.15f);
			backPlane.material = grid;
			
			// A sphere for testing
			CSGSphere sphere = new CSGSphere();
			sphere.material = refractive;
			
			// Collect objects in intersectable list
			IntersectableList intersectableList = new IntersectableList();

			
			intersectableList.add(sphere);
			intersectableList.add(backPlane);
			
			// Set the root node for the scene
			root = intersectableList;
			
			// Light sources
			Vector3f lightPos = new Vector3f(eye);
			lightPos.add(new Vector3f(-1.f, 0.f, 0.f));
			LightGeometry pointLight1 = new PointLight(lightPos, new Spectrum(14.f, 14.f, 14.f));
			lightPos.add(new Vector3f(2.f, 0.f, 0.f));
			LightGeometry pointLight2 = new PointLight(lightPos, new Spectrum(14.f, 14.f, 14.f));
			LightGeometry pointLight3 = new PointLight(new Vector3f(0.f, 7.f, 0.f), new Spectrum(14.f, 14.f, 14.f));
			lightList = new LightList();
			lightList.add(pointLight1);
			lightList.add(pointLight2);
			lightList.add(pointLight3);
		}
		
	}
	<beautiful scenes>+=
	new RefractiveSphere(),

		
<h2>Appendix: Utilities</h2>
The dot product can be used to determine whether two vectors n and d lie in the same halfspace:
<img src=dothalfspace.JPG></img>
	<determine halfspace>=
	public static boolean sameHalfspace(Vector3f n, Vector3f d) {
		return n.dot(d) > 0;
	}
	
The following function computes ray-sphere intersection times 
(<a href=http://sci.tuomastonteri.fi/programming/sse/example3>source</a>). 
Note: we can return all hit points,
also the ones with negative t-value, that is, points that lie "behind"
the origin of the ray.
Use dFac -1 to solve for first, 1 for second hitpoint.
    <ray sphere intersection>=
		public static float intersectSphere(Vector3f center, float radius, Ray r, float dFac) {
			float a = r.direction.dot(r.direction); // Squared length of ray direction vector. > 0.
			float b = r.direction.dot(M.scale(2.0f, M.sub(r.origin, center)));
			float c = center.dot(center) + r.origin.dot(r.origin) 
					- 2.0f*r.origin.dot(center)
					- radius*radius;
			// ^^ end of vector stuff
			
			// Discriminant of quadratic equation.
			float D = b*b + (-4.0f)*a*c;
			
			// If ray can not intersect then stop
			if (D < 0) return Float.NEGATIVE_INFINITY;
			D = (float)Math.sqrt(D);

			// Ray can intersect the sphere, solve the closer hitpoint (negative D)
			float a2 = 2 * a;
			float t = -b/a2 + dFac*D/a2;
			return t;
		}
	[rt/M.java]= 
	package rt;
	<common imports>
	public class M {
        public static final float PI = (float)Math.PI;
		<determine halfspace>
		/** Squared distance between v1 and v2, ||v1 - v2||^2 */
		public static float dist2(Tuple3f v1, Tuple3f v2)
		{
			return sub(v1, v2).lengthSquared();
		}
		
		/** v1 - v2, Vector that points from v2 to v1. */
		public static Vector3f sub(Tuple3f v1, Tuple3f v2)
		{
			Vector3f r = new Vector3f(v1);
			r.sub(v2);
			return r;
		}
		
		public static Vector3f normalize(Vector3f v) {
			Vector3f r = new Vector3f(v);
			r.normalize();
			return r;
		}
		
		public static float clamp(float f, float min, float max) {
			return f < min ? min : f > max ? max : f; // TODO would Math.min/max be faster?
		}

		public static float clamp(float f) {
			return clamp(f, 0.f, 1.f);
		}
		
		public static float sqrtf(float f) {
			return (float)Math.sqrt(f);
		}
		
        public static float sqr(float f) {
            return f*f;
		}

		public static float expf(float a) {
			return (float)Math.exp(a);
		}

		public static float powf(float a, float b) {
			return (float)Math.pow(a, b);
		}
		
		/** v1 + v2 */
		public static Vector3f add(Tuple3f v1, Tuple3f v2)
		{
			Vector3f r = new Vector3f(v1);
			r.add(v2);
			return r;
		}
		
		/** v1 x v2, vector that is perpendicular to both v1 and v2 and has length = area of parallelogramm
		 * spanned by the two. */
		public static Vector3f cross(Vector3f v1, Vector3f v2)
		{
			Vector3f r = new Vector3f();
			r.cross(v1, v2);
			return r;
		}
		
		/** - v */
		public static Vector3f negate(Vector3f v)
		{
			Vector3f r = new Vector3f(v);
			r.negate();
			return r;
		}
		
		/** s (Scalar) * v */
		public static Vector3f scale(float s, Vector3f v)
		{
			Vector3f r = new Vector3f(v);
			r.scale(s);
			return r;
		}
        
        public static Vector3f scale(Vector3f v, float s) {return scale(v,s);}
		
		/** m^-1 */
		public static Matrix4f invert(Matrix4f m)
		{
			Matrix4f r = new Matrix4f(m);
			r.invert();
			return r;
		}

		public static Vector3f transformVectorAsPoint(Matrix4f t, Vector3f v) {
			Point3f p = new Point3f(v);
			t.transform(p);
			return new Vector3f(p);
		}

		/** o + d * t */
		public static Vector3f t(Vector3f o, Vector3f d, float t) {
			Vector3f n = new Vector3f(d);
			n.scaleAdd(t, o);
			return n;
		}

		<math utilities>
	}
		
A utility class to help us run benchmarks.
	[rt/Timer.java]= 
	package rt;
    <common imports>
	public class Timer {
		private long _startTime;
        
		public Timer() { _startTime = this.timeNow(); }
		
		public long timeElapsed() {
			return this.timeNow() - _startTime;
		}

		private long timeNow() {return new Date().getTime();}
	} 	
	
<h2>Appendix: Obj Reader</h2>
The only external scene and data description format that we currently support are obj files.
This reads an .obj file including normals and stores it in a Mesh.
scale scales the object to fit into a cube of the given size
	[rt/ObjReader.java]= 
	package rt;
	<common imports>
	public class ObjReader {
		public static Mesh read(String fileName, float scale) {
			Mesh mesh;
			try {
				mesh = ObjReader._read(fileName, scale);
			} catch(IOException e) {
				throw new RuntimeException("Could not read .obj file: "+fileName);
			}
			return mesh;
		}
		
		private static Mesh _read(String fileName, float scale) throws IOException
		{
			BufferedReader reader;
			ArrayList<float[]> vertices = new ArrayList<float[]>();
			ArrayList<float[]> texCoords = new ArrayList<float[]>();
			ArrayList<float[]> normals = new ArrayList<float[]>();
			ArrayList<int[][]> faces = new ArrayList<int[][]>();
			
			boolean hasNormals, hasTexCoords;
			hasNormals = true;
			hasTexCoords = true;
			
			// Extents for normalization
			float xMin, xMax, yMin, yMax, zMin, zMax;
			xMin = Float.MAX_VALUE;
			xMax = Float.MIN_VALUE;
			yMin = Float.MAX_VALUE;
			yMax = Float.MIN_VALUE;
			zMin = Float.MAX_VALUE;
			zMax = Float.MIN_VALUE;
			
			reader = new BufferedReader(new FileReader(fileName));

			String line = null;
			while((line = reader.readLine()) != null)
			{	
				// Read line
				String[] s = line.split("\\s+");
				
				// Parse
				if(s[0].compareTo("v")==0)
				{
					// Position
					float[] v = new float[3];
					v[0] = Float.valueOf(s[1]).floatValue();
					v[1] = Float.valueOf(s[2]).floatValue();
					v[2] = Float.valueOf(s[3]).floatValue();
					vertices.add(v);
					
					// Update extent
					if(v[0] < xMin) xMin = v[0];
					if(v[0] > xMax) xMax = v[0];
					if(v[1] < yMin) yMin = v[1];
					if(v[1] > yMax) yMax = v[1];
					if(v[2] < zMin) zMin = v[2];
					if(v[2] > zMax) zMax = v[2];
				} 
				else if(s[0].compareTo("vn")==0)
				{
					// Normal
					float[] n = new float[3];
					n[0] = Float.valueOf(s[1]).floatValue();
					n[1] = Float.valueOf(s[2]).floatValue();
					n[2] = Float.valueOf(s[3]).floatValue();
					normals.add(n);
				}
				else if(s[0].compareTo("vt")==0)
				{
					// Texture
					float[] t = new float[2];
					t[0] = Float.valueOf(s[1]).floatValue();
					t[1] = Float.valueOf(s[2]).floatValue();
					texCoords.add(t);
				}
				else if(s[0].compareTo("f")==0)
				{
					// Indices
					int[][] indices = new int[3][3];
					
					// For all vertices
					int i=1;
					while(i < s.length)
					{	
						// Get indices for vertex position, tex. coords., and normals
						String[] ss = s[i].split("/");
						int k=0;
						while(k < ss.length)
						{
							if(ss[k].length()>0)
								indices[i-1][k] = Integer.valueOf(ss[k]).intValue();
							else
							{
								indices[i-1][k] = -1;
								if(k == 1) hasTexCoords = false;
								if(k == 2) hasNormals = false;
							}
							k++;
						}
						if(ss.length == 1)
						{
							hasTexCoords = false;
							hasNormals = false;
						}
						i++;
					}
					faces.add(indices);
				}
				else if(s[0].length()>0 && s[0].charAt(0)!='#')
				{
                    // usemtl statements etc generate this
					//System.err.print("Unknown token '".concat(line).concat("'\n"));
				}
			}

			// Normalization
			float xTrans = -(xMax+xMin)/2;
			float yTrans = -(yMax+yMin)/2;
			float zTrans = -(zMax+zMin)/2;
			float xScale = 2/(xMax-xMin);
			float yScale = 2/(yMax-yMin);
			float zScale = 2/(zMax-zMin);
			float s = yScale;
			if(xScale < yScale) s = xScale;
			if(zScale < s) s = zScale;
			scale = s*scale;
			
			// Brute force approach to generate single index per vertex
			// Expand arrays
			int nFaces = faces.size();
			float[] verticesFinal = new float[nFaces*9];
			float[] normalsFinal = new float[nFaces*9];
			float[] texCoordsFinal = new float[nFaces*6];
			int[] indices = new int[nFaces*3];
			
			// For all faces
			int vertexNr = 0;
			for(int i=0; i<nFaces; i++)
			{
				// For all vertices
				for(int j=0; j<3; j++)
				{
					// Copy positions, tex. coords., and normals to expanded arrays
					// Note: we subtract one from the index because indexing in the obj
					// file is 1-based, whereas our arrays are 0-based
					verticesFinal[vertexNr*3] = vertices.get(faces.get(i)[j][0]-1)[0];
					verticesFinal[vertexNr*3+1] = vertices.get(faces.get(i)[j][0]-1)[1];
					verticesFinal[vertexNr*3+2] = vertices.get(faces.get(i)[j][0]-1)[2];
					
					verticesFinal[vertexNr*3] = scale*(verticesFinal[vertexNr*3]+xTrans);
					verticesFinal[vertexNr*3+1] = scale*(verticesFinal[vertexNr*3+1]+yTrans);
					verticesFinal[vertexNr*3+2] = scale*(verticesFinal[vertexNr*3+2]+zTrans);
					
					if(hasNormals)
					{
						normalsFinal[vertexNr*3] = normals.get(faces.get(i)[j][2]-1)[0];
						normalsFinal[vertexNr*3+1] = normals.get(faces.get(i)[j][2]-1)[1];
						normalsFinal[vertexNr*3+2] = normals.get(faces.get(i)[j][2]-1)[2];
					} 
					
					if(hasTexCoords)
					{
						texCoordsFinal[vertexNr*2] = texCoords.get(faces.get(i)[j][1]-1)[0];
						texCoordsFinal[vertexNr*2+1] = texCoords.get(faces.get(i)[j][1]-1)[1];
					}
					
					indices[vertexNr] = vertexNr;
					vertexNr++;
				}
				
				if(!hasNormals)
				{
					Vector3f d0 = new Vector3f(
						verticesFinal[(vertexNr-1)*3  ]-verticesFinal[(vertexNr-3)*3],
						verticesFinal[(vertexNr-1)*3+1]-verticesFinal[(vertexNr-3)*3+1],
						verticesFinal[(vertexNr-1)*3+2]-verticesFinal[(vertexNr-3)*3+2]);
					Vector3f d1 = new Vector3f(
						verticesFinal[(vertexNr-2)*3  ]-verticesFinal[(vertexNr-3)*3],
						verticesFinal[(vertexNr-2)*3+1]-verticesFinal[(vertexNr-3)*3+1],
						verticesFinal[(vertexNr-2)*3+2]-verticesFinal[(vertexNr-3)*3+2]);
					Vector3f n = new Vector3f();
					n.cross(d1,d0);
					n.normalize();
					for(int j=0; j<3; j++)
					{
						normalsFinal[(vertexNr-(j+1))*3]   = n.x;
						normalsFinal[(vertexNr-(j+1))*3+1] = n.y;
						normalsFinal[(vertexNr-(j+1))*3+2] = n.z;
					}
				}
			}

			reader.close();
			return new Mesh(verticesFinal, normalsFinal, indices);
		}
	}