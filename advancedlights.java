

<h5>Physically based shading</h5>
In the study of such algorithms, the following concepts are of central concern:
<h6>Surface Tangents</h6> 
The normalized tangent vectors at the hit point are 
computed from the normal on-demand, when constructTangents() is called.
t1, t2, normal is a right handed frame.
	<hit record datastructure>+=
		public Vector3f t1, t2;
		
		public void constructTangents()  {
			t1 = M.tangentTo(normal);
			t2 = M.cross(normal, t1);
		}

<p>
Materials implement functionality for shading surfaces using their BRDFs.

Light sources are implemented using materials that return an emission term.
	<material methods and data>+=
		public Spectrum evaluateEmission(HitRecord hitRecord, Vector3f wOut);
evaluates emission for outgoing direction. This method is typically called 
by an integrator when the integrator obtained the outgoing direction of
the emission by sampling a point on a light source.

@param hitRecord Information about hit point on light source
@param wOut Outgoing direction, normalized and pointing away from the surface


		
		/**
		 * Evaluate BRDF for pair of incoming and outgoing directions. This method
		 * is typically called by an integrator when the integrator obtained the incident 
		 * direction by sampling a point on a light source.
		 * 
		 * @param hitRecord Information about hit point
		 * @param wOut Outgoing direction, normalized and pointing away from the surface
		 * @param wIn Incoming direction, normalized and pointing away from the surface
		 * @return BRDF value
		 */
		public Spectrum evaluateBRDF(HitRecord hitRecord, Vector3f wOut, Vector3f wIn);

        
		/**
		 * Return whether material has perfect specular reflection. 
		 */
		public boolean hasSpecularReflection();
		
        <material methods and data>+=
        <shading sample>
		/**
		 * Evaluate specular reflection. This method is typically called by a recursive
		 * ray tracer to follow the path of specular reflection.
		 */
		public ShadingSample evaluateSpecularReflection(HitRecord hitRecord);
		
		/**
		 * Return whether the material has perfect specular refraction. 
		 */
		public boolean hasSpecularRefraction();

		/**
		 * Evaluate specular refraction. This method is typically called by a recursive
		 * ray tracer to follow the path of specular refraction.
		 */
		public ShadingSample evaluateSpecularRefraction(HitRecord hitRecord);	
		
		/**
		 * Calculate a shading sample, given a uniform random sample as input. This 
		 * method is typically called in a path tracer to sample and evaluate the
		 * next path segment. The methods decides which component of the material to 
		 * sample (diffuse, specular reflection or refraction, etc.), computes an 
		 * incident direction, and returns the BRDF value, the direction, and the 
		 * probability density (stored in a {@link ShadingSample}). 
		 */
		public ShadingSample getShadingSample(HitRecord hitRecord, float[] sample);

		/**
		 * Calculate an emission sample, given a hit record and a uniform random 
		 * sample as input. This method is typically called in a bidirectional
		 * path tracer to sample and evaluate the first light path segment. The 
		 * methods computes an outgoing direction, and returns the emission value, 
		 * the direction, and the probability density (all stored in a 
		 * {@link ShadingSample}). 
		 */
		public ShadingSample getEmissionSample(HitRecord hitRecord, float[] sample);

		/**
		 * Indicates whether the material casts shadows or not. 
		 */
		public boolean castsShadows();
	}

		

<h3>Materials with procedural colors</h3>	

	
<h3>Shading Sample</h3>
A "shading sample" is a brdf (bidirectional reflectance distribution function) value, an emission value,
a sampled direction, information on wether this is a specular sample and the probability density of 
the sample.
	<shading sample>=
	public class ShadingSample {
		public Spectrum brdf;
		public Spectrum emission;
		
		/**
		 * The sampled direction.
		 */
		public Vector3f w;
		
		/**
		 * Tells the integrator whether this is a specular sample. In this case,
		 * a cosine factor in the specular BRDF should be omitted in the returned 
		 * BRDF value, and the integrator should act accordingly.
		 */
		public boolean isSpecular;
		
		/**
		 * The (directional) probability density of the sample
		 */
		public float p;
		
		public ShadingSample(Spectrum brdf, Spectrum emission, Vector3f w, boolean isSpecular, float p)
		{
			this.brdf = new Spectrum(brdf);
			this.emission = new Spectrum(emission);
			this.w = new Vector3f(w);
			this.isSpecular = isSpecular;
			this.p = p;
		}
		
		public ShadingSample()
		{			
		}
	}	
    


<h4>Point Light Shadowless Integrator</h4>
Given a two dimensional random uniformly distributed sample in (0,1)^2,
this should return a uniformly sampled location on the surface of the intersectable.
    <additional intersectable methods and data>+=
        public Vector3f sample(float[] s) {
            throw new RuntimeException(this.getClass() + " does not implement sample() method");
        }

This is easy to implement for a rectangle:
    <additional intersectable methods and data>+=
        public Vector3f sample(float[] s) {
            return M.add(a, M.add(M.scale(s[0], da), M.scale(s[1], db)));
        }
        
Basic morphological shading integrator with diffuse and specular terms
computed for any amount of lights. It simply iterates over the light sources and accumulates
their contributions. 
No shadow testing, reflection, refraction, or area light sources, etc. is supported.

	<intersect scene>=
        HitRecord hitRecord = root.intersect(r);
        if (hitRecord == null) return new Spectrum();
        
We immediately return a „background color“ if nothing was hit.
	<compute contribution s of lightsource>=
        <sample lightsource>
        <get direction and distance to lightsource>
        <evaluate brdf for pointlight>
        
	<evaluate brdf for pointlight>=
        Spectrum s = hitRecord.material.evaluateBRDF(hitRecord, hitRecord.w, lightDir); 
        s.mult(<color of light source>);
				
We assume we are dealing with point lights only, and we let the integrator take care of the quadratic falloff of the strength of such lightsources with distance.
	<evaluate brdf for pointlight>+=					
        s.mult(1.f/d2);
				
	<sample lightsource>=
        HitRecord lightHit = lightSource.sample(make2dSample());
        
    <get direction and distance to lightsource>=
        Vector3f lightDir = M.sub(lightHit.position, hitRecord.position);
        float d2 = lightDir.lengthSquared();
        lightDir.normalize();
        
	<color of light source>=
        lightHit.material.evaluateEmission(lightHit, M.negate(lightDir))		

	[rt/integrators/PointLightIntegrator.java]= 
	package rt.integrators;

	<common imports>

	public class PointLightIntegrator extends Integrator {

		LightList lightList;
		Intersectable root;
		
		public PointLightIntegrator(Scene scene)
		{
            super(scene);
			this.lightList = scene.getLightList();
			this.root = scene.getIntersectable();
		}

		public Spectrum integrate(Ray r) {
			<intersect scene>

			Spectrum outgoing = new Spectrum(0.f, 0.f, 0.f);	
			for (LightGeometry lightSource : lightList) {	
				<compute contribution s of lightsource>
				outgoing.add(s);				
			}			
			return outgoing;		
		}

	}
		
	[rt/integrators/PointLightIntegratorFactory.java]= 
	package rt.integrators;
	<common imports>
	public class PointLightIntegratorFactory extends IntegratorFactory {
		public Integrator make(Scene scene) {return new PointLightIntegrator(scene);}
	}
		
<h4>Whitted Integrator</h4>
Extension of the point light integrator that supports shadows, and recursive raytracing to implement reflection and refraction.
Since the latter two features might result in an infinite loop, we limit the recursion depth.
	<recursion depth limit>=
	5

	<evaluate specular reflection and refraction>=
	if (depth < <recursion depth limit> 
		// && M.sameHalfspace(hitRecord.w, hitRecord.normal) // Only reflect/refract on front surfaces
		) {
		depth++;
		Spectrum ref = new Spectrum();
		Spectrum refl = null, refr = null;
		
		if (hitRecord.material.hasSpecularReflection())
			refl = integrate(hitRecord, hitRecord.material.evaluateSpecularReflection(hitRecord), depth);
			
		if (hitRecord.material.hasSpecularRefraction())
			refr = integrate(hitRecord, hitRecord.material.evaluateSpecularRefraction(hitRecord), depth);

		<combine reflection and refraction>

		outgoing.add(ref);
	} 

<img src=fresnel.png></img>			
	<combine reflection and refraction>=
	if (refl != null) {
		ref.add(refl);
		if (refr != null) {
			float F = ((Refractive)hitRecord.material).schlickF(hitRecord);
			refl.mult(F);
			refr.mult(1-F);
			ref = refl.add(refr);
		}
	} else if (refr != null)
		ref.add(refr);

As for the lighting, we just need to add the shadowtest, adding some bias to avoid shadow acne
	<shadow test>=
	HitRecord shadowHit = root.intersect(
		new Ray(M.t(hitRecord.position, lightDir, BIAS), lightDir));
	if (shadowHit != null && 
		<hit before lightsource> &&
		shadowHit.material.castsShadows()) 
		<no contribution from this lightsource>

	<hit before lightsource>=
        shadowHit.t*shadowHit.t <= d2
    
	<no contribution from this lightsource>=
        continue;
		
	[rt/integrators/WhittedIntegrator.java]= 
	package rt.integrators;
	<common imports>
	public class WhittedIntegrator extends PointLightIntegrator {

		public WhittedIntegrator(Scene scene){super(scene);}
		
		/** Only ShadingSample’s direction w is considered (is that how shadingSample should be used?) */
		public Spectrum integrate(HitRecord r, ShadingSample s, int depth) {
			if (s == null) return null;
			Vector3f p = new Vector3f(s.w); p.scale(BIAS); p.add(r.position);
			return integrate(new Ray(p, s.w), depth+1);
		}
		
		public Spectrum integrate(Ray r) {
			return integrate(r, 0);
		}
		
		public Spectrum integrate(Ray r, int depth) {
			<intersect scene>

			Spectrum outgoing = new Spectrum();
			<evaluate specular reflection and refraction>
			for (LightGeometry lightSource : lightList) {	
				<sample lightsource>
				<get direction and distance to lightsource>
				<shadow test>
				<evaluate brdf for pointlight>
				
				outgoing.add(s);
			}
			
			return outgoing;	
		}
	}
	
	[rt/integrators/WhittedIntegratorFactory.java]= 
	package rt.integrators;
	<common imports>
	public class WhittedIntegratorFactory extends IntegratorFactory {
		public Integrator make(Scene scene) {return new WhittedIntegrator(scene);}
	}
    
<h4>Area Light Integrator</h4>
Extension of the point light integrator that supports soft shadows
by randomly sampling an area light source
    <light visibility rays>=
	64
times. 	
	[rt/integrators/AreaLightIntegrator.java]= 
	package rt.integrators;
	<common imports>
	public class AreaLightIntegrator extends WhittedIntegrator {
        public AreaLightIntegrator(Scene s) {super(s);}
		public Spectrum integrate(Ray r, int depth) {
			<intersect scene>

			Spectrum outgoing = new Spectrum();
			<evaluate specular reflection and refraction>
			for (LightGeometry lightSource : lightList) {	
                for (int i = 0; i < <light visibility rays>; i++) {
                    <sample lightsource>
                    <get direction and distance to lightsource>
                    <shadow test>
                    <evaluate brdf for pointlight>
                    
                    s.mult(1.f/<light visibility rays>);
                    outgoing.add(s);
                }
			}
			
			return outgoing;	
		}
	}
	
	[rt/integrators/AreaLightIntegratorFactory.java]= 
	package rt.integrators;
	<common imports>
	public class AreaLightIntegratorFactory extends IntegratorFactory {
		public Integrator make(Scene scene) {return new AreaLightIntegrator(scene);}
	}
    
<h2>Lights</h2>
An interface to implement light sources. Light sources derive from 
this interface, and they store a reference to a Material
with an emission term.
The method
	<sample light point>+=
	public abstract HitRecord sample(float[] s);
should return a random point on the surface of the light source given a random value s in [0..1].
	[rt/LightGeometry.java]= 
	package rt;
	<common imports>
	public abstract class LightGeometry extends Intersectable {
		<sample light point>
	}
	[rt/LightList.java]= 
	package rt;
	import java.util.ArrayList;
	public class LightList extends ArrayList<LightGeometry> {}
		
<h3>PointLight</h3>
Implements a point light using a PointLightMaterial.
	[rt/lightsources/PointLight.java]= 
	package rt.lightsources;

	<common imports>
	public class PointLight extends LightGeometry {

		Vector3f position;
		PointLightMaterial pointLightMaterial;
		Random rand;
		
		public PointLight(Vector3f position, Spectrum emission)
		{
			this.position = new Vector3f(position);
			this.rand = new Random();
			pointLightMaterial = new PointLightMaterial(emission);
		}
		public HitRecord intersect(Ray r) {return null;}

		public HitRecord sample(float[] s) {
			HitRecord hitRecord = new HitRecord();
			hitRecord.position = new Vector3f(position);
			hitRecord.material = pointLightMaterial;
			hitRecord.p = 1.f;
			return hitRecord;
		}
	}
	
	[rt/materials/PointLightMaterial.java]= 
	package rt.materials;

	<common imports>
	public class PointLightMaterial extends Material {

		Spectrum emission;
		Random rand;
		
		public PointLightMaterial(Spectrum emission)
		{
			this.emission = new Spectrum(emission);
			this.rand = new Random();
		}
		
		public Spectrum evaluateEmission(HitRecord hitRecord, Vector3f wOut) {
			return new Spectrum(emission);
		}

		/**
		 * Return a random direction over the full sphere of directions.
		 */
		public ShadingSample getEmissionSample(HitRecord hitRecord, float[] sample) {
			// To be implemented
			return null;
		}
		
		// The rest of the functions should not be called
		public ShadingSample getShadingSample(HitRecord hitRecord, float[] sample) {<should not be called>}
		public boolean castsShadows() {
			<should not be called>}
		public Spectrum evaluateBRDF(HitRecord hitRecord, Vector3f wOut,
				Vector3f wIn) {<should not be called>}
		public boolean hasSpecularReflection() {<should not be called>}
		public ShadingSample evaluateSpecularReflection(HitRecord hitRecord) {<should not be called>}
		public boolean hasSpecularRefraction() {<should not be called>}
		public ShadingSample evaluateSpecularRefraction(HitRecord hitRecord) {<should not be called>}
	}
	<should not be called>=
	throw new RuntimeException();

<h3>AreaLight</h3>
Implements an area light using an AreaLightMaterial.
	[rt/lightsources/AreaLight.java]= 
	package rt.lightsources;

	<common imports>
	public class AreaLight extends LightGeometry {
        Rectangle r; // TODO this should be a LightGeometry
		Material m;
		
		public AreaLight(Vector3f a, Vector3f da, Vector3f db, Spectrum emission)
		{
			this.r = new Rectangle(a,da,db);
			m = new AreaLightMaterial(emission);
		}
		public HitRecord intersect(Ray ray) {
            HitRecord h = r.intersect(ray);
            if (h!=null) h.intersectable = this;
            return h;
        }

		public HitRecord sample(float[] s) {
			HitRecord hitRecord = new HitRecord();
			hitRecord.position = r.sample(s);
			hitRecord.material = m;
			hitRecord.p = 1.f;
			return hitRecord;
		}
	}
	
	[rt/materials/AreaLightMaterial.java]= 
	package rt.materials;

	<common imports>
	public class AreaLightMaterial extends Material {

		Spectrum emission;
		Random rand;
		
		public AreaLightMaterial(Spectrum emission)
		{
			this.emission = new Spectrum(emission);
			this.rand = new Random();
		}
		
		public Spectrum evaluateEmission(HitRecord hitRecord, Vector3f wOut) {
			return new Spectrum(emission);
		}

		/**
		 * Return a random direction over the full sphere of directions.
		 */
		public ShadingSample getEmissionSample(HitRecord hitRecord, float[] sample) {
			// To be implemented
			return null;
		}
		
		// The rest of the functions should not be called
		public ShadingSample getShadingSample(HitRecord hitRecord, float[] sample) {<should not be called>}
		public boolean castsShadows() {
			<should not be called>}
		public Spectrum evaluateBRDF(HitRecord hitRecord, Vector3f wOut,
				Vector3f wIn) {<should not be called>}
		public boolean hasSpecularReflection() {<should not be called>}
		public ShadingSample evaluateSpecularReflection(HitRecord hitRecord) {<should not be called>}
		public boolean hasSpecularRefraction() {<should not be called>}
		public ShadingSample evaluateSpecularRefraction(HitRecord hitRecord) {<should not be called>}
	}
	

Same thing but with the uv debug integrator.
	[rt/testscenes/ObjectUVTest.java]= 
	package rt.testscenes;
	<common imports>
	public abstract class ObjectUVTest extends ObjectTest {
        public ObjectUVTest() {this(new Vector3f());}
		public ObjectUVTest(Vector3f eye) {
			super(eye);
			integratorFactory = new UVDebugIntegratorFactory();
		}
	}