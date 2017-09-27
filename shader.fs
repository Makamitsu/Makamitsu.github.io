#define M_PI 3.1415926535897932384626433832795
precision mediump float;

//Uniform for the cam.
uniform vec3 cam_origin;
uniform float cam_focal;
uniform float cam_width;
uniform float cam_height;
uniform float cam_rotationPhi;
uniform float cam_rotationTheta;

uniform vec3 cam_direction;
uniform vec2 screenSize;
//
//Uniform for LIGHT_1.
uniform vec3 light1_center;
uniform vec3 light1_power;
//
//Uniform for LIGHT_2.
uniform vec3 light2_center;
uniform vec3 light2_power;
//
//Uniform for LIGHT_3.
uniform vec3 light3_center;
uniform vec3 light3_power;
//
//Uniform for SPHERE_1.
uniform vec3 sphere1_center;
uniform float sphere1_radius;
uniform vec3 sphere1_kd;
uniform float sphere1_ks;
uniform float sphere1_m;
uniform float sphere1_ni;
//
//Uniform for SPHERE_2.
uniform vec3 sphere2_center;
uniform float sphere2_radius;
uniform vec3 sphere2_kd;
uniform float sphere2_ks;
uniform float sphere2_m;
uniform float sphere2_ni;
//
//Uniform for SPHERE_3.
uniform vec3 sphere3_center;
uniform float sphere3_radius;
uniform vec3 sphere3_kd;
uniform float sphere3_ks;
uniform float sphere3_m;
uniform float sphere3_ni;
//
//Uniform for SPHERE_4.
uniform vec3 sphere4_center;
uniform float sphere4_radius;
uniform vec3 sphere4_kd;
uniform float sphere4_ks;
uniform float sphere4_m;
uniform float sphere4_ni;
//
//Uniform for SPHERE_5.
uniform vec3 sphere5_center;
uniform float sphere5_radius;
uniform vec3 sphere5_kd;
uniform float sphere5_ks;
uniform float sphere5_m;
uniform float sphere5_ni;
//
//Uniform for PLAN_1.
uniform vec3 plan1_normal;
uniform float plan1_shift;
uniform vec3 plan1_kd;
uniform float plan1_ks;
uniform float plan1_m;
uniform float plan1_ni;
//
//Uniform for PLAN_2.
uniform vec3 plan2_normal;
uniform float plan2_shift;
uniform vec3 plan2_kd;
uniform float plan2_ks;
uniform float plan2_m;
uniform float plan2_ni;
//
//Uniform for PLAN_3.
uniform vec3 plan3_normal;
uniform float plan3_shift;
uniform vec3 plan3_kd;
uniform float plan3_ks;
uniform float plan3_m;
uniform float plan3_ni;
//
//Uniform for PLAN_4.
uniform vec3 plan4_normal;
uniform float plan4_shift;
uniform vec3 plan4_kd;
uniform float plan4_ks;
uniform float plan4_m;
uniform float plan4_ni;
//
//Uniform for PLAN_5.
uniform vec3 plan5_normal;
uniform float plan5_shift;
uniform vec3 plan5_kd;
uniform float plan5_ks;
uniform float plan5_m;
uniform float plan5_ni;
//


const int nb_light = 3;
const int nb_sphere = 5;
const int nb_plan = 5;

//DEFINITION des STRUCT
struct Ray {
    vec3 o;
    vec3 v;
};
//
struct Material {
    vec3 kd;
    float ks;
    float m;
    float ni;
    float damier;
    float mirror;
};
//
struct Light {
    vec3 c;
    vec3 p;
};
//
struct Sphere {
    vec3 c;
    float r;
    Material mat;
};
//
struct Plan {
    vec3 n;
    float s;
    Material mat;
};
struct Intersect {
    float t;
    vec3 n;
    Material mat;
};

Light lights[nb_light];
Sphere spheres[nb_sphere];
Plan plans[nb_plan];

void createScene(){
    Light light1;
    light1.c = light1_center;
    light1.p = light1_power;

    Light light2;
    light2.c = light2_center;
    light2.p = light2_power;

    Light light3;
    light3.c = light3_center;
    light3.p = light3_power;

    Sphere sphere1;
    sphere1.c = sphere1_center;
    sphere1.r = sphere1_radius;
    Material s1mat;
    s1mat.kd = sphere1_kd;
    s1mat.ks = sphere1_ks;
    s1mat.m = sphere1_m;
    s1mat.ni = sphere1_ni;
    sphere1.mat = s1mat;

    Sphere sphere2;
    sphere2.c = sphere2_center;
    sphere2.r = sphere2_radius;
    Material s2mat;
    s2mat.kd = sphere2_kd;
    s2mat.ks = sphere2_ks;
    s2mat.m = sphere2_m;
    s2mat.ni = sphere2_ni;
    sphere2.mat = s2mat;

    Sphere sphere3;
    sphere3.c = sphere3_center;
    sphere3.r = sphere3_radius;
    Material s3mat;
    s3mat.kd = sphere3_kd;
    s3mat.ks = sphere3_ks;
    s3mat.m = sphere3_m;
    s3mat.ni = sphere3_ni;
    sphere3.mat = s3mat;

    Sphere sphere4;
    sphere4.c = sphere4_center;
    sphere4.r = sphere4_radius;
    Material s4mat;
    s4mat.kd = sphere4_kd;
    s4mat.ks = sphere4_ks;
    s4mat.m = sphere4_m;
    s4mat.ni = sphere4_ni;
    sphere4.mat = s4mat;

    Sphere sphere5;
    sphere5.c = sphere5_center;
    sphere5.r = sphere5_radius;
    Material s5mat;
    s5mat.kd = sphere5_kd;
    s5mat.ks = sphere5_ks;
    s5mat.m = sphere5_m;
    s5mat.ni = sphere5_ni;
    s5mat.mirror = 0.2;
    sphere5.mat = s5mat;

    Plan plan1;
    plan1.n = plan1_normal;
    plan1.s = plan1_shift;
    Material p1mat;
    p1mat.kd = plan1_kd;
    p1mat.ks = plan1_ks;
    p1mat.m = plan1_m;
    p1mat.ni = plan1_ni;
    plan1.mat = p1mat;

    Plan plan2;
    plan2.n = plan2_normal;
    plan2.s = plan2_shift;
    Material p2mat;
    p2mat.kd = plan2_kd;
    p2mat.ks = plan2_ks;
    p2mat.m = plan2_m;
    p2mat.ni = plan2_ni;
    plan2.mat = p2mat;

    Plan plan3;
    plan3.n = plan3_normal;
    plan3.s = plan3_shift;
    Material p3mat;
    p3mat.kd = plan3_kd;
    p3mat.ks = plan3_ks;
    p3mat.m = plan3_m;
    p3mat.ni = plan3_ni;
    plan3.mat = p3mat;

    Plan plan4;
    plan4.n = plan4_normal;
    plan4.s = plan4_shift;
    Material p4mat;
    p4mat.kd = plan4_kd;
    p4mat.ks = plan4_ks;
    p4mat.m = plan4_m;
    p4mat.ni = plan4_ni;
    plan4.mat = p4mat;

    Plan plan5;
    plan5.n = plan5_normal;
    plan5.s = plan5_shift;
    Material p5mat;
    p5mat.kd = plan5_kd;
    p5mat.ks = plan5_ks;
    p5mat.m = plan5_m;
    p5mat.ni = plan5_ni;
    p5mat.mirror = 0.05;
    plan5.mat = p5mat;

    lights[0] = light1;
    lights[1] = light2;
    lights[2] = light3;
    spheres[0] = sphere1;
    spheres[1] = sphere2;
    spheres[2] = sphere3;
    spheres[3] = sphere4;
    spheres[4] = sphere5;
    plans[0] = plan1;
    plans[1] = plan2;
    plans[2] = plan3;
    plans[3] = plan4;
    plans[4] = plan5;
}

//functions de calcul matriciel
mat4 lookat(vec3 eye, vec3 target, vec3 up)
{
    vec3 zaxis = normalize(eye - target);
    vec3 xaxis = normalize(cross(up, zaxis));
    vec3 yaxis = cross(zaxis, xaxis);

    mat4 orientation = mat4(
       xaxis[0], yaxis[0], zaxis[0], 0,
       xaxis[1], yaxis[1], zaxis[1], 0,
       xaxis[2], yaxis[2], zaxis[2], 0,
         0,       0,       0,     1);

    mat4 translation = mat4(
              1,       0,       0, 0,
              0,       1,       0, 0,
              0,       0,       1, 0,
        -eye[0], -eye[1], -eye[2], 1);

    return orientation;// * translation;
}

vec3 rotationX (vec3 v,float a){
    return v * mat3(cos(a), -sin(a), 0.0,
             sin(a), cos(a), 0.0,
             0.0,0.0,1.0);
}
vec3 rotationY (vec3 v,float a){
    return v * mat3(cos(a),0.0, -sin(a),
             0.0,1.0,0.0,
             sin(a), 0.0, cos(a));
}
vec3 rotationZ (vec3 v,float a){
    return v * mat3(1.0, 0.0, 0.0,
             0.0,cos(a), -sin(a),
             0.0,sin(a), cos(a));
}
mat3 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}
vec3 rotationPhi(vec3 v, float a){
    return v * mat3(cos(a), 0,-sin(a),
        0, 1,0,
        sin(a),0,cos(a));
}

vec3 rotationTheta(vec3 v, float a){
    return v * mat3(cos(a),sin(a),0,
        -sin(a), cos(a), 0,
        0, 0, 1 );
}

mat4 lookAtMatrix(vec3 origin, vec3 target){
    vec3 forward;
    vec3 rigth;
    vec3 up;
    mat4 matrix2;
    mat4 resutlMatrix;

    forward = normalize(origin - target);
    rigth = normalize(cross(forward, vec3(0.0,1.0,0.0)));
    up = cross(rigth, forward);

    matrix2 = mat4(rigth.x , up.x, forward.x, 0.0,
                   rigth.y , up.y, forward.y, 0.0,
                   rigth.z , up.z, forward.z, 0.0,
                   0.0, 0.0, 0.0, 1.0);

   return matrix2;
}


float intersectionSphere(Sphere s, Ray r) {
	vec3 o = r.o - s.c;
	float a = dot(r.v,r.v);
	float b = 2.0*dot(r.v,o);
	float c = dot(o,o) - s.r*s.r;
	float delta = b*b -(4.0*a*c);


    if(delta < 0.0){
        return -1.0;
    }else{
        if(delta == 0.0){
            return -b/2.0*a;
        }else{
            float t1 = (-b-sqrt(delta))/(2.0*a);
            float t2 = (-b+sqrt(delta))/(2.0*a);
            if(t1 < 0.0){
                return t2;
            }
            if(t2 < 0.0){
                return t1;
            }
            return min(t1,t2);
        }
    }
}

float intersectionPlan(Plan p , Ray r){
     return -(dot(p.n,r.o)+p.s)/(dot(p.n,r.v));
}

Intersect throw(Ray r) {
    float tAux;
    Intersect inter;
    inter.t = -1.0;

    for(int i = 0; i < nb_sphere; i++){
        tAux = intersectionSphere(spheres[i], r);
        if (tAux > 0.0 && tAux < inter.t || inter.t == -1.0 && tAux != -1.0 && tAux > 0.0){
            inter.t = tAux;
            inter.mat = spheres[i].mat;
            inter.n = normalize(((r.v * tAux) + r.o) - spheres[i].c); //VECTEUR NORMAL (impact - centre sphere)
        }
    }
    for(int i = 0; i < nb_plan; i++){
        tAux = intersectionPlan(plans[i], r);
        if (tAux > 0.0 && tAux < inter.t || inter.t == -1.0 && tAux != -1.0 && tAux > 0.0){
            inter.t = tAux;
            inter.mat = plans[i].mat;
            inter.n = plans[i].n; //VECTEUR NORMAL le meme que le plan
        }
    }
    return inter;
}

float getG (vec3 n, vec3 Vi, vec3 V0){
	vec3 h = normalize (Vi+V0);
	float G = min(1.0,min(2.0*dot(n,h)*dot(n,V0)/dot(V0,h),2.0*dot(n,h)*dot(n,Vi)/dot(V0,h) ) );
	return (G);
}
float getDistribution (vec3 n, vec3 Vi, vec3 V0, Material mat){
	vec3 h = normalize (Vi+V0);
	float cosalpha = dot(n,h);
	float D = (1.0 / mat.m*pow(cosalpha,4.0)) * exp(-pow((tan(acos(cosalpha))/mat.m),2.0));
	return (D);
	}
float getFresnel (vec3 n, vec3 Vi, vec3 V0, Material mat){
	vec3 h = normalize (Vi+V0);
	float c = dot(V0,h);
	float g = sqrt( pow(mat.ni,2.0) + pow(c,2.0) - 1.0 );
	float F = 0.5*(pow((g-c),2.0))/(pow((g+c),2.0))*(1.0+((pow(c*(g+c)-1.0,2.0))/(pow(c*(g-c)+1.0,2.0))));
	return (F);
}
vec3 microFacette(Light Li, vec3 n, vec3 Vi, vec3 V0, Material mat){
	float costeta =  dot(Vi,n);
	float D  = getDistribution ( n,  Vi,  V0,  mat);
	float F = getFresnel ( n,  Vi,  V0,  mat);
	float G = getG ( n,  Vi,  V0);
	return vec3( Li.p * ((mat.kd/M_PI) + (D*F*G / 4.0 * dot(Vi,n)*dot(V0,n)) ) * costeta );
}

vec3 blendLights(Intersect impact, Ray r){
    if(impact.t > 0.0){ //t doit etre positif et !=1.0

        Intersect isObscuration;

        vec3 i = r.o + r.v*impact.t;

        vec3 BRDF = vec3(0.0);

        for(int l = 0; l < nb_light; l++){
            vec3 vi = lights[l].c-i;//vecteur vers la lumière

            Ray toLight;
            toLight.v = normalize(vi);
            toLight.o = i + normalize(vi)/10.0;

            isObscuration = throw(toLight);

            if(isObscuration.t < 0.0 ||     // Si le rayon vers la lumière ne touche aucun élément (< 0.0 ou == -1.0)
            length(isObscuration.t*toLight.v) > length(lights[l].c-toLight.o)){// OU si l'element touché est derriere la lumière
                vec3 v0 = normalize(r.o-i); //Vecteur vers l'origine
                vi = normalize(vi);

                if(impact.mat.damier > 0.0){
                    if( mod(i.z,2.0) < impact.mat.damier){
                        impact.mat.kd*0.0;
                        BRDF += microFacette(lights[l], impact.n , vi, v0, impact.mat);
                    }else{
                        BRDF += microFacette(lights[l], impact.n , vi, v0, impact.mat);
                    }

                }else{
                    BRDF += microFacette(lights[l], impact.n , vi, v0, impact.mat);
//                    BRDF -= vec3(0.05);
                }

            }
        }


    return BRDF;

    }else {
        return vec3(0.8); //Couleur du fond.
    }
}

void main(void){
//creation de la scene
    createScene();

//creation du rayon
    Ray ray;
	ray.v =
	normalize(
	vec3(
        cam_width*(gl_FragCoord.x-screenSize.x/2.0)/(screenSize.x/2.0),
        cam_height*(gl_FragCoord.y-screenSize.y/2.0)/(screenSize.y/2.0),
        cam_focal
        )
    );
    ray.o = cam_origin;

//application des matrices de rotations

    mat4 viewMatrix = lookat(ray.o, ray.o+vec3(cam_direction.x, 0.0, cam_direction.z), vec3(0.0,1.0,0.0));
    ray.v = (viewMatrix * vec4(ray.v,1.0) ).xyz;


//lancement du rayon et retour d'un élément intersect au point d'impact
    Intersect impact = throw(ray);


    vec3 BRDF = blendLights(impact, ray);

    vec3 BRDF2;

    if (impact.mat.mirror > 0.0){
        Ray ray2;
        ray2.v = reflect(ray.v,impact.n);
        ray2.o = (ray.o+impact.t*ray.v)+ray2.v/100.0;
        Intersect impact2 = throw(ray2);

        BRDF2 = blendLights(impact2, ray2);
    }else{
        BRDF2 = BRDF;
    }



//calcule de la lumiere locale au point d'impact

// VALIDATION DE LA NORMALE AU POINT D'IMPACT
//    gl_FragColor = vec4(impact.n,1.0);
//    gl_FragColor = vec4(impact.n/2.0+0.5,1.0);
// VALIDATION DU TMIN DE L'ORIGINE AU POINT D'IMPACT
//    gl_FragColor = vec4(vec3(impact.t)/100.0,1.0);
// VALIDATION DE LA RACUPERATION DES MATERIAUX
//    gl_FragColor = vec4(impact.mat.kd,1.0);

    gl_FragColor = vec4(((1.0-impact.mat.mirror)*BRDF+impact.mat.mirror*BRDF2)/2.0,1.0);

}

/**
* Ce qui a ete verifier
* fonction throw()
*       calculs des normales (plans + spheres + ombres) = ok!
*
*
*/