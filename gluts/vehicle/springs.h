#ifndef _ENG_PHYSICS_SPRINGS__H__
#define _ENG_PHYSICS_SPRINGS__H__

#ifndef EXTERN_INLINE
#define EXTERN_INLINE static inline
#endif

/*
    Integrates a spring that conforms to
    
	force = -(k*m) * x - (b*m) * velocity

    or

	accel = - k*x - b*velocity

    input is current position and velocity
*/
EXTERN_INLINE float physicsSpringGetVelocity(float x, float velocity, float dt, float k, float b);


/*
    Integrate a spring that conforms to
    
	velocity(t) = -spring * x(t)

    Which has an expontential solution x(t) = Amplitude * expf(-spring * x(t) )
    NB: this spring is equivalient to lerping
*/
EXTERN_INLINE float physicsExpSpringGetVelocity(float x, float spring, float dt);


/*
==================================================================================  
    INLINES
==================================================================================
*/

EXTERN_INLINE float physicsSpringGetVelocity(float x, float velocity, float dt, float k, float b)
{
    float a = dt*k;
    return (velocity - a*x) / (1.f + dt*(a + b));
}


EXTERN_INLINE float physicsExpSpringGetVelocity(float x, float spring, float dt)
{
#if 0
    return (x/dt) * (expf(-spring * dt) - 1.f);
#else
    // Cheaper implicit solution version
    return (-x * spring)/(1.f+spring*dt);
#endif
}


#endif
