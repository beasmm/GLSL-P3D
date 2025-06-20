/**
* ver hash functions em
* https://www.shadertoy.com/view/XlGcRh hash functions GPU
* http://www.jcgt.org/published/0009/03/02/
 */

 #include "./common.glsl"
 #iChannel0 "self"
 
#define SCENE 0
#define POINT false
#define QUAD true
#define MICROFACETS false

bool hit_world(Ray r, float tmin, float tmax, inout HitRecord rec)
{
    bool hit = false;
    rec.t = tmax;

    #if SCENE == 0       //Shirley Weekend scene

        if(hit_quad(createQuad(vec3(-10.0, -0.05, 10.0), vec3(10.0, -0.05, 10.0), vec3(10.0, -0.05, -10.0), vec3(-10.0, -0.05, -10.0)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2));
        }

        if(hit_sphere(createSphere(vec3(-4.0, 1.0, 0.0), 1.0), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2, 0.95, 0.1));
            //rec.material = createDiffuseMaterial(vec3(0.4, 0.2, 0.1));
        }

        if(hit_sphere(createSphere(vec3(4.0, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            if (MICROFACETS)
                rec.material = createPlasticMaterial(vec3(0.2), vec3(0.7, 0.6, 0.5), 0.04, 0.0);
            else 
                rec.material = createMetalMaterial(vec3(0.7, 0.6, 0.5), 0.0);
        }

        if(hit_sphere(createSphere(vec3(-1.5, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0), 1.33, 0.0);
        }

        if(hit_sphere(createSphere(vec3(-1.5, 1.0, 0.0), -0.5),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0), 1.33, 0.0);
        }

        if(hit_sphere(createSphere(vec3(1.5, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0, 0.9, 0.9), 1.5, 0.0);
        }
            
        int numxy = 5;
        
        for(int x = -numxy; x < numxy; ++x)
        {
            for(int y = -numxy; y < numxy; ++y)
            {
                float fx = float(x);
                float fy = float(y);
                float seed = fx + fy / 1000.0;
                vec3 rand1 = hash3(seed);
                vec3 center = vec3(fx + 0.9 * rand1.x, 0.2, fy + 0.9 * rand1.y);
                float chooseMaterial = rand1.z;
                if(distance(center, vec3(4.0, 0.2, 0.0)) > 0.9)
                {
                    if(chooseMaterial < 0.3)
                    {
                        vec3 center1 = center + vec3(0.0, hash1(gSeed) * 0.5, 0.0);
                        // diffuse
                        if(hit_movingSphere(createMovingSphere(center, center1, 0.2, 0.0, 1.0),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.5)
                    {
                        // diffuse
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.7)
                    {
                        // metal
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, 0.0);
                        }
                    }
                    else if(chooseMaterial < 0.9)
                    {
                        // metal
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, hash1(seed));
                        }
                    }
                    else
                    {
                        // glass (Dielectric)
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDielectricMaterial(hash3(seed), 1.33, 0.0);
                        }
                    }
                }
            }
        }
    #elif SCENE == 1 //from https://blog.demofox.org/2020/06/14/casual-shadertoy-path-tracing-3-fresnel-rough-refraction-absorption-orbit-camera/

        // diffuse floor
        
            vec3 A = vec3(-25.0f, -12.5f, 10.0f);
            vec3 B = vec3( 25.0f, -12.5f, 10.0f);
            vec3 C = vec3( 25.0f, -12.5f, -5.0f);
            vec3 D = vec3(-25.0f, -12.5f, -5.0f);

            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7));
            }

        //stripped background
        {
            vec3 A = vec3(-25.0f, -10.5f, -5.0f);
            vec3 B = vec3( 25.0f, -10.5f, -5.0f);
            vec3 C = vec3( 25.0f, -1.5f, -5.0f);
            vec3 D = vec3(-25.0f, -1.5f, -5.0f);
        
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                float shade = floor(mod(rec.pos.x, 1.0f) * 2.0f);
                rec.material = createDiffuseMaterial(vec3(shade));
            }
        }

        // ceiling piece above light
        
        {
            vec3 A = vec3(-7.5f, 12.5f, 5.0f);
            vec3 B = vec3( 7.5f, 12.5f, 5.0f);
            vec3 C = vec3( 7.5f, 12.5f, -5.0f);
            vec3 D = vec3(-7.5f, 12.5f, -5.0f);

            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7));
            }
        }    
       
        // light
        
        {
            vec3 A = vec3(-5.0f, 12.3f,  2.5f);
            vec3 B = vec3( 5.0f, 12.3f,  2.5f);
            vec3 C = vec3( 5.0f, 12.3f,  -2.5f);
            vec3 D = vec3(-5.0f, 12.3f,  -2.5f);

             if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.0));
                rec.material.emissive = vec3(1.0f, 0.9f, 0.9f) * 20.0f;
            }
        }
 
        const int c_numSpheres = 7;
        for (int sphereIndex = 0; sphereIndex < c_numSpheres; ++sphereIndex)
        {
            vec3 center = vec3(-18.0 + 6.0 * float(sphereIndex), -8.0, 0.0);
            if(hit_sphere(createSphere(center, 2.8),r,tmin,rec.t,rec))
            {
                hit = true;
                float r = float(sphereIndex) / float(c_numSpheres-1) * 0.1f;
                rec.material = createDielectricMaterial(vec3(0.0, 0.5, 1.0), 1.1, r);
            }
        }

    #elif SCENE == 2
    #elif SCENE == 3
    #endif

    return hit;
}

vec3 directlighting(pointLight pl, Ray r, HitRecord rec) {
    vec3 toLight = normalize(pl.pos - rec.pos);
    vec3 normal = normalize(rec.normal);
    vec3 viewDir = normalize(-r.d); // View direction is the inverse of the ray

    // Shadow ray
    Ray shadowRay;
    shadowRay.o = rec.pos + 0.001 * normal;
    shadowRay.d = toLight;

    HitRecord dummy;
    if (hit_world(shadowRay, 0.001, length(pl.pos - rec.pos), dummy)) {
        // In shadow
        return vec3(0.0);
    }

    // Diffuse component
    float diff = max(dot(normal, toLight), 0.0);
    vec3 diffCol = diff * rec.material.albedo * pl.color;

    // Specular component
    vec3 reflectDir = reflect(-toLight, normal);
    float shininess = 32.0; // Use a constant shininess value
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), shininess);
    vec3 specColor = vec3(1.0); // Use white as the specular color
    vec3 specCol = spec * specColor * pl.color;

    return diffCol + specCol;
}

// new code for quad light
vec3 directLightingFromQuadLight(QuadLight ql, Ray r, HitRecord rec) {
    // Sample a random point on the quad light
    vec3 lightPos = ql.quad.a + hash3(gSeed) * (ql.quad.b - ql.quad.a) + 
                   hash3(gSeed) * (ql.quad.d - ql.quad.a);
    
    vec3 toLight = normalize(lightPos - rec.pos);
    float distance = length(lightPos - rec.pos);
    
    // Shadow ray
    Ray shadowRay;
    shadowRay.o = rec.pos + 0.001 * rec.normal;
    shadowRay.d = toLight;

    HitRecord dummy;
    if (hit_world(shadowRay, 0.001, distance, dummy)) {
        return vec3(0.0); // In shadow
    }

    // Calculate the diffuse component
    float diff = max(dot(rec.normal, toLight), 0.0);
    
    // Calculate light area for proper normalization
    vec3 edge1 = ql.quad.b - ql.quad.a;
    vec3 edge2 = ql.quad.d - ql.quad.a;
    float area = length(cross(edge1, edge2));
    
    // Solid angle approximation
    float lightCos = max(dot(-toLight, normalize(cross(edge1, edge2))), 0.0);
    float solidAngle = (lightCos * area) / (distance * distance);
    
    return diff * ql.color * ql.intensity * solidAngle;
}


#define MAX_BOUNCES 10

vec3 rayColor(Ray r)
{
    HitRecord rec;
    vec3 col = vec3(0.0);
    vec3 throughput = vec3(1.0);

    for(int i = 0; i < MAX_BOUNCES; ++i)
    {
        if(hit_world(r, 0.001, 10000.0, rec))
        {
            // === 1. Add direct lighting from three point lights ===
            if (POINT)
            {
                pointLight pl0 = createPointLight(vec3(-10.0, 15.0, 0.0), vec3(1.0));
                pointLight pl1 = createPointLight(vec3(8.0, 15.0, 3.0), vec3(1.0));
                pointLight pl2 = createPointLight(vec3(1.0, 15.0, -9.0), vec3(1.0));

                col += throughput * directlighting(pl0, r, rec);
                col += throughput * directlighting(pl1, r, rec);
                col += throughput * directlighting(pl2, r, rec);
            }

            if (QUAD) 
            {
                Quad q = Quad(
                    vec3(-2.0, 5.0, -2.0),
                    vec3(2.0, 5.0, -2.0),   
                    vec3(2.0, 5.0, 2.0),   
                    vec3(-2.0, 5.0, 2.0)
                );
                QuadLight ql = createQuadLight(q, vec3(1.0, 0.9, 0.8), 1.0);

                col += throughput * directLightingFromQuadLight(ql, r, rec);
            }

            col += throughput * rec.material.emissive;

            Ray scatterRay;
            vec3 atten;
            if(scatter(r, rec, atten, scatterRay))
            {
                throughput *= atten;
                r = scatterRay;
            }
            else
            {
                break;
            }
        }
        else
        {
            float t = 0.8 * (r.d.y + 1.0);
            vec3 skyColor = mix(vec3(1.0), vec3(0.5, 0.7, 1.0), t);
            col += throughput * skyColor;
            break;
        }
    }

    return col;
}


#define MAX_SAMPLES 10000.0

void main()
{
    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + iTime;

    vec2 mouse = iMouse.xy / iResolution.xy;
    mouse.x = mouse.x * 2.0 - 1.0;

    vec3 camTarget = vec3(0.0);

    // make sure camera is above the base plane
    float minElevation = 0.05; // 3° above horizontal
    float maxElevation = 0.6;  // 34°
    float elevation = mix(minElevation, maxElevation, mouse.y);

    float azimuth = 2.0 * 3.14159 * mouse.x; // horizontal angle

    // Clamp zoom radius to a reasonable range
    float radius = mix(5.0, 15.0, 1.0 - mouse.y);

    vec3 camPos = camTarget + radius * vec3(
        cos(elevation) * sin(azimuth),
        sin(elevation),
        cos(elevation) * cos(azimuth)
    );

    float fovy = 60.0;
    float aperture = 8.0;
    float distToFocus = 1.0;
    float time0 = 0.0;
    float time1 = 1.0;
    Camera cam = createCamera(
        camPos,
        camTarget,
        vec3(0.0, 1.0, 0.0),    // world up vector
        fovy,
        iResolution.x / iResolution.y,
        aperture,
        distToFocus,
        time0,
        time1);

//usa-se o 4 canal de cor para guardar o numero de samples e não o iFrame pois quando se mexe o rato faz-se reset

    vec4 prev = texture(iChannel0, gl_FragCoord.xy / iResolution.xy);
    vec3 prevLinear = toLinear(prev.xyz);  

    vec2 ps = gl_FragCoord.xy + hash2(gSeed);
    //vec2 ps = gl_FragCoord.xy;
    vec3 color = rayColor(getRay(cam, ps));

    if(iMouseButton.x != 0.0 || iMouseButton.y != 0.0)
    {
        gl_FragColor = vec4(toGamma(color), 1.0);  //samples number reset = 1
        return;
    }
    if(prev.w > MAX_SAMPLES)   
    {
        gl_FragColor = prev;
        return;
    }

    float w = prev.w + 1.0;
    color = mix(prevLinear, color, 1.0/w);
    gl_FragColor = vec4(toGamma(color), w);
}
