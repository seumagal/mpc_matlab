#define V          xC[0]
#define U          xC[1]
#define YAW_RATE   xC[2]
#define ROLL_RATE  xC[3]
#define ROLL       xC[4]

#define FL (0)
#define FR (1)
#define RL (2)
#define RR (3)

#define MFA	 (&params[0])
#define MFB  (&params[8])
#define Z    (&params[17]) 
#define ZF0  (params[24]) 
#define ZR0  (params[25])
#define A    (params[26])
#define B    (params[27])
#define TF2   (params[28])
#define TR2   (params[29])
#define D    (&params[30])

#define LAT_ACCEL  ( dstates[0] )
#define LONG_ACCEL ( dstates[1] )

#define sind(x) (sin(fmod((x),360) * M_PI / 180))
#define tand(x) (tan(fmod((x),360) * M_PI / 180))
#define atand(x) (fmod(atan(x)*180/M_PI,360) )
#define atan2d(x,y) (fmod(atan2(x,y)*180/M_PI,360) )

real_T slip_angle[4] ={0,0,0,0};
real_T long_f[4] ={0,0,0,0};
real_T lat_f[4] ={0,0,0,0};
real_T load[4] ={0,0,0,0};
real_T dstates[4] = {0,0,0,0};
const real_T long_slip[4] = {0.1,0.1,0.1,0.1};
