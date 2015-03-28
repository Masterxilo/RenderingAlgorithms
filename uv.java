

This shows the hitpoint uv values.
<img src=uvs.png></img>
	[rt/integrators/UVDebugIntegrator.java]= 
	package rt.integrators;
	<common imports>
	public class UVDebugIntegrator extends DebugIntegrator {
		public UVDebugIntegrator(Scene scene) {super(scene);}
		
		public Spectrum integrate(Ray r) {
			HitRecord hitRecord = scene.getIntersectable().intersect(r);
			if (hitRecord == null) return new Spectrum(0.f,0.f,0.f);
			if (hitRecord.t <= 0.f) return new Spectrum(1.f,0.f,0.f);
			return new Spectrum(hitRecord.u, hitRecord.v, 0);
		}

	}
		
	[rt/integrators/UVDebugIntegratorFactory.java]= 
	package rt.integrators;
	import rt.*;
	public class UVDebugIntegratorFactory extends IntegratorFactory {
		public Integrator make(Scene scene) {return new UVDebugIntegrator(scene);}
	}