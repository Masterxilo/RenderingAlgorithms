
<h4>Light Simulation</h4>
Interesting and more (photo) realistic images are obtained once we add 
a true simulation of light and light-transport to the scene.

<h5>Distinction of light sources</h5>
Some of the techniques illustrated below require that we 
can distinguish between objects that 
emit light into the scene from those which do not.
We also need to be able to enumerate all of these.
We can achieve this by simply adding a flag to intersectables and collecting 
those intersectables with the flag set into a list of "light sources".
    <additional intersectable methods and data>+=
    public boolean isLight = false;
    
    <additional scene data>+=
		protected IntersectableList lightList;
		public IntersectableList getLightList() {
            return lightList;
        }

<h3>Incident Ray/To-Viewer/Eye/Observer/Out Direction w</h3>
The direction towards the origin
of the ray that hit the surface is often useful for shading. 
This unit vector points away from the surface, in the direction
opposite to that of the incident ray.
	<hit record datastructure>+=
		public Vector3f w() {return M.negate(ray.direction);}

<h3>Surface Tangents</h3> 
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

		
<h3>Blinn-Phong Shading Model</h3>
Note that the parameter value kd is the diffuse reflectance,
which should be in the range [0,1], a value of 1 meaning all light
is reflected (diffusely), and none is absorbed.
The diffuse BRDF corresponding to kd is actually kd/pi.
	<normalize kd>=
	this.kd.mult(1/(float)Math.PI); 

ks is the specular color, which is white for most real materials.
s is the shininess parameter.
<img src=shade.png></img>
The Blinn BRDF formula computes the reflected light for a given direction towards the light and towards the observer.
L is also called the incident direction.
	<blinn evaluate brdf>=
	public Spectrum evaluateBRDF(HitRecord hitRecord, Vector3f e, Vector3f L) {
		Vector3f n = hitRecord.normal;
		Vector3f h = <unit halfway vector>;
		return 
			new Spectrum(kd).mult(M.clamp(L.dot(n))).add(
			new Spectrum(ks).mult(M.powf (h.dot(n), s))
		);
	}

<img src=blinn.png></img>
	<unit halfway vector>=
	M.normalize(M.add(L, e))

	[rt/materials/Blinn.java]= 
	package rt.materials;
	<common imports>

	public class Blinn extends Material {
		float s;
		Spectrum ks;
		Spectrum kd;
		public Blinn(Spectrum kd, Spectrum ks, float s) {
			this.kd = kd; this.ks = ks; this.s = s;
			<normalize kd>
		}

		<blinn evaluate brdf>

		public boolean castsShadows() {return true;}

		public boolean hasSpecularReflection() {return false;}
		public ShadingSample evaluateSpecularReflection(HitRecord hitRecord) {return null;}
		public boolean hasSpecularRefraction() {return false;}
		public ShadingSample evaluateSpecularRefraction(HitRecord hitRecord) {return null;}
		public ShadingSample getShadingSample(HitRecord hitRecord, float[] sample) {return null;}
		public Spectrum evaluateEmission(HitRecord hitRecord, Vector3f wOut) {return new Spectrum(0.f, 0.f, 0.f);}
		public ShadingSample getEmissionSample(HitRecord hitRecord, float[] sample) {return new ShadingSample();}
	}

A basic diffuse material has no (black) specular color.
	[rt/materials/Diffuse.java]= 
	package rt.materials;
	<common imports>
	
	public class Diffuse extends Blinn {
		public Diffuse(Spectrum kd) {
			super(kd, new Spectrum(), 0);
		}

		public Diffuse() {this(new Spectrum(1.f, 1.f, 1.f));}
	}

<h3>Reflective Material</h3>
A direction d reflected on the normal can be constructed as:
	<hit record datastructure>+=
	public Vector3f reflect(Vector3f d) {
		return M.sub(M.scale(2 * d.dot(normal), normal), d);
	}
		
Most of the time, we will need to reflect w, the incident direction, e.g. to implement the law of reflection.
	<hit record methods>+=
	public Vector3f reflectedW() { return reflect(w()); }
    
This material is used to represent mirrors.
The reflection is easy to do (we only use the w component of ShadingSample).
	<define reflection>=
	public boolean hasSpecularReflection() {return true;}
	public ShadingSample evaluateSpecularReflection(HitRecord hitRecord) {
		Vector3f d = hitRecord.reflectedW();
		return new ShadingSample(new Spectrum(), new Spectrum(), d, true, 1.f);
	}
	
	[rt/materials/Reflective.java]= 
	package rt.materials;
	<common imports>
	public class Reflective extends Material 
	{
		<define reflection>
			
		public boolean hasSpecularRefraction() {return false;}
		public ShadingSample evaluateSpecularRefraction(HitRecord hitRecord) {return null;}
	
		public Spectrum evaluateBRDF(HitRecord hitRecord, Vector3f wOut,
				Vector3f wIn) {return new Spectrum();}
		public Spectrum evaluateEmission(HitRecord hitRecord, Vector3f wOut) {return new Spectrum();}
		public ShadingSample getShadingSample(HitRecord hitRecord, float[] sample) {return null;}
		public ShadingSample getEmissionSample(HitRecord hitRecord, float[] sample) {return null;}
		public boolean castsShadows() {return true;}
	}
	
<h3>Refractive Material</h3>
<h3>The 'Hit Plane'</h3>	
For constructing rays that lie within the plane 
spanned by the normal and the incident direction 
w (call this the hit plane), it is useful to have the surface's
tangent that lies in this plane. 
We call this vector ipt and compute it on-demand.
	<hit record datastructure>+=
        private Vector3f ipt;
        private void updateIpt() {
            assert normal != null : "normal must be non-null";
            assert w != null : "w must be non-null";
            Vector3f cnw = M.cross(normal, w); 
            ipt = cnw; ipt.cross(normal, cnw);
            ipt.normalize();
        }
		
ipt and the normal then form a cartesian coordinate system on the hit plane, on which we can now construct arbitrary points.
	<hit record datastructure>+=
		public Vector3f hitPlanePoint(float x, float y) {
			updateIpt();
			Vector3f r = M.scale(x, ipt);
			r.add(M.scale(y, normal));
			return r;
		}
		
We can of course give these in polar coordinates.
	<hit record datastructure>+=
		public Vector3f hitPlanePointPolar(float ang, float len) {
			return hitPlanePoint((float)Math.cos(ang)*len, (float)Math.sin(ang)*len);
		}
		public Vector3f hitPlanePointPolar(float ang) {
			return hitPlanePointPolar(ang, 1.f);
		}
		
This will be used for refraction.
In the following image, the black arrow is the surface normal, red is cnw, green is ipt, purple is w. 
In blue, we show a vector constructed using polar coordinates on the ipt-normal plane (hit plane).
<img src=refex.png></img>

    <unit tests>+=
    @Test
    public void testHitPlanePointPolar() {
        HitRecord h = new HitRecord(
            new Ray(new Vector3f(), new Vector3f(1.f, -1.f, 0)), 
            0, null
            new Vector3f(0.f, 1.f, 0)
        ));
        Vector3f o = h.hitPlanePointPolar(M.PI/4.f);
        assertEquals(0.707107f, o.x, 0.0001f);
        assertEquals(0.707107f, o.y, 0.0001f);
        assertEquals(0, o.z, 0.0001f);
    }
    
    @Test
    public void testHitPlanePointPolar2() {
        HitRecord h = new HitRecord(
            new Ray(new Vector3f(), new Vector3f(1.f, 1.f, 0)), 
            0, null
            new Vector3f(0.f, 1.f, 0)
        ));
        Vector3f o = h.hitPlanePointPolar(M.PI/4.f);
        assertEquals(0.707107f, o.x, 0.0001f);
        assertEquals(0.707107f, o.y, 0.0001f);
        assertEquals(0, o.z, 0.0001f);
    }

To determine the refracted angle, 
we need to be able to tell the incident angle.
	<hit record datastructure>+=
	public float angleToNormal(Vector3f v) {
		v = M.normalize(v);
        assert normal != null : "normal must be non-null";
		double a = (double)v.dot(normal);
		return (float)Math.acos(a);
	}
Computes the angle between the given vector and the normal of the surface, 
assuming both vectors point to the same hemisphere.
v need not be normalized.
    <unit tests>+=
    @Test
    public void testAngleToNormal() {
        HitRecord h = new HitRecord(null, 0, null, new Vector3f(1.f, 0, 0));
        float w = h.angleToNormal(new Vector3f(1.f, 1.f, 0));
        assertEquals(Math.PI/4.f, w, 0.0001f);
    }

<h4>The Refractive Material</h4>
This material is used to represent glass, water and other refractive/reflective materials.



In physics, the „amount“ of refraction of a material is characterized by its refractive index n. 
As a parameter for our surface material here, n is the refractive index of the medium below the surface (pointing in the direction of the normal) using this material.
Refraction is computed using the incoming angle of the light and the refractive indices on both sides of the surface according to Snell’s law:
<img src=snell.png></img>

	<define refraction>=
	public boolean hasSpecularRefraction() {return true;}
	public ShadingSample evaluateSpecularRefraction(HitRecord hitRecord) {
		<determine incoming angle>
		<determine medium and adjust angle>
		
        assert th1 >= 0 && th1 <= M.PI/2;
		float a = (float)Math.sin(th1)*n1/n2;
        
		<check total internal reflection>
		float th2 = (float)Math.asin(a);
		assert th2 >= 0 && th2 <= M.PI/2;
		<compute thout>

		return new ShadingSample(
            new Spectrum(), 
            new Spectrum(), 
            hitRecord.hitPlanePointPolar(thout), true, 1.f);
	}
    
We only use the ‚direction‘ value of the ShadingSample (I have no idea what the rest is for).
We make this method return null if there should be total internal reflection.
	<check total internal reflection>=
	if (a > 1) return null; 
We also have
	<determine incoming angle>=
	float th1 = hitRecord.angleToNormal(hitRecord.w);
and
	<determine medium and adjust angle>=
	float n1, n2;
	boolean leaving = false;
	if (th1 > Math.PI/2) { // leaving medium 
		n1 = n; n2 = 1.f;
		th1 = (float)Math.PI - th1;
        leaving = true;
	} else {
		n1 = 1.f; n2 = n;
	}

determines whether we are leaving or entering the medium represented by this material’s n and sets n1 and n2 accordingly. The „outside“ medium is assumed to be air (vacuum actually).

<p>
The „outgoing angle“ thout is measured, by definition of the .hitPlanePointPolar method, relative to the +x axis in the ipt, n coordinate system in the hit-plane. This is the line surrounded by n1, v1, n2, v2 in the above picture. Since th2 is measured against the normal (or against -normal if we are entering the medium) we need to adjust accordingly 
	<compute thout>=
	float thout = leaving ? (M.PI/2 - th2) : (-M.PI/2 + th2);
	
In the following image, the black arrow is the surface normal, green is ipt (green, black form the x and y axis of the coordinate system against which thout is measured), red is cnw, purple is w (see definition of hitPlanePointPolar). 
Blue is the refracted vector.
<img src=refex.png></img>
If w comes from the direction opposite to the normal, the angle of the refracted vector must be measured against the non-inverted normal.
<img src=refex2.png></img>

<img src=fresnel.png></img>
Gives the interpolation factor in [0,1] to be used between the 
reflected and refracted color.
	<schlick fresnel factor>=
	public float schlickF(HitRecord hitRecord) {
		<determine incoming angle>
		<determine medium and adjust angle>
			
		float f = (1-n1/n2)*(1-n1/n2) / ((1+n1/n2)*(1+n1/n2));
		return f + (1-f)*M.powf(1-hitRecord.w.dot(hitRecord.normal), 5);
	}

<h4>Overall</h4>
	[rt/materials/Refractive.java]= 
	package rt.materials;
	<common imports>
	public class Refractive extends Reflective
	{
		float n;
		public Refractive(float n) {this.n = n;}

		<define refraction>
		<schlick fresnel factor>
	}
	
We can also define a material that, for artistic purposes, has reflection only (such materials do not exist).
	[rt/materials/RefractiveOnly.java]= 
	package rt.materials;
	<common imports>
	public class RefractiveOnly extends Refractive {
		public RefractiveOnly(float n) {super(n);}
		public boolean hasSpecularReflection() {return false;}
        public ShadingSample evaluateSpecularReflection(HitRecord hitRecord) {return null;}
	}
	
Let us refract a ray coming from -1,-1,0 on the yz plane at a surface with refractive index 2.
    <unit tests>+=
    @Test
    public void testRefract() {
        HitRecord h = new HitRecord();
        h.w = new Vector3f(1.f, 1.f, 0); // to camera
        h.normal = new Vector3f(1.f, 0.f, 0);
        Material m = new RefractiveOnly(2.f);
        
        Vector3f o = m.evaluateSpecularRefraction(h).w;
        
        assertEquals(-0.935414f, o.x, 0.0001f);
        assertEquals(-0.353553f, o.y, 0.0001f);
        assertEquals(o.z, 0, 0.0001f);
    }
    
    @Test
    public void testRefract2() {
        HitRecord h = new HitRecord();
        h.w = new Vector3f(-1.f, 1.f, 0); // to camera
        h.normal = new Vector3f(1.f, 0.f, 0);
        Material m = new RefractiveOnly(1.3f);
        
        Vector3f o = m.evaluateSpecularRefraction(h).w;
         
        assertEquals(0.3937f, o.x, 0.0001f);
        assertEquals(-0.919239f, o.y, 0.0001f);
        assertEquals(o.z, 0, 0.0001f);
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
		static final float BIAS = 0.0001f;
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