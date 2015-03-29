
	
<h3>Beautiful Scenes</h3>


    
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
 
	

