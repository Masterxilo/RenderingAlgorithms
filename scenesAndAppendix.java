
	
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

