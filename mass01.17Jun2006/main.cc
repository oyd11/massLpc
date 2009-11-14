/*
 * Mass-Spring system
 * Attached to (prob buggy) Kelly-Lochbaum lattice filter
 *
 * For Peter Pabon's Physical Modeling class, May 2006
 *
 * By JK Calev (Oyd11)
 * http://oyd11.planet-d.net/
 *  oyd111
 *  on the gmail server
 *
 *
 *  Uses SDL (SimpleDirectmediaLayer) for IO
 *  'should' compiler under win32/linux/macOSX
 *  (tested on win32 only)
 *
 *  @:+ - marks 'highlights' in the code (search for it)
 */

#include <stdlib.h>

#include <math.h>

#include "SDL/SDL.h"

#define PI M_PI
#define PI_2 (PI*2)


#define SZ_X 700
#define SZ_Y 420

#define min(X,Y) (X)>(Y)?(Y):(X)
// de'f in some cpp header as a templ as well

struct Mass {
    double x,v,a,k;
};

#define MAX_MASS 64
Mass tract[MAX_MASS];
// TODO: could be more convenient to use 3 'double' arrays
// and manage 'x,v,a' values 'manually'

//int nMass = 54;
int nMass = 51;
double invSR=1./5000.; // springSystem SR
// TODO: impl KR4 - and reduce sampling rate - should work!
int mouseColumn=-1;

double k=0.4;
double damp=1.; // damping... in [0..1]
bool skipMass[MAX_MASS];

double sc=100.; // scale from simulation-space to pixels
double minSize=0.4;
double limit_size=true;
double plusSize;

double tGlobal=0.; // just some counter...

// AUDIO vars:
//int SR = 11025; // SamplingRate
int SR = 11025*4; // SamplingRate
double SF; //= 1./SR; // can this be trusted in compiletime?
int aBuffSz=2048; // buffSize in bytes or samples?
double tractLen = 17.; // cm
double cmPerSample; // speedOfSound/SR

// interface vars:
bool pause=false;
SDL_Surface *screen;

void updateTract();

void init_masses(void){
    for (int i=0;i<MAX_MASS;++i) {
	tract[i].v=tract[i].a=tract[i].x=0.;
	skipMass[i]=false;
	tract[i].k=0.9;
    }
    tract[1].x=-0.9;     
    tract[2].x=-0.3;     
    tract[3].x=-0.1;     
    tract[7].x=-0.1;     
    tract[8].x=-0.2;     

    tract[19].x=-0.1;     
    tract[23].x=-0.3;     
    tract[29].x=-0.2;     

    skipMass[1]=skipMass[nMass-2]=true;
    tract[22].x=-0.3;     
    skipMass[22]=true;
    for (int i=23;i<33;++i) {
	tract[i].k=0.4;
    }

    tract[32].x=-0.32;     
    skipMass[32]=true;
}

double mousePull[MAX_MASS];

// @:+  MassSpring system:
void game_tick(void) {
    // 'leaky integrator'
    for (int i= 1 ; i<nMass-1; ++i ) {
	if (skipMass[i]) continue;
	//      t=t+dt; y=y+v*dt; v=v+a*dt; a=-k*y/m; //??
	tract[i].v = tract[i].v*damp+invSR*tract[i].a;
	tract[i].x += invSR*tract[i].v;
	tract[i].a =-(tract[i].x-tract[i+1].x)*tract[i+1].k
	    +(tract[i-1].x-tract[i].x)*tract[i].k;
//	-(tract[i].x-1.4*mousePull[i]);
	// XXX: bug!!??
	// pull abit more than 'mousePull', for "mouse" (or 'vowel')
	// to be the target
    }
}


void reset(){
}

void putpixel(unsigned int x,unsigned int y, int color) {
    if (x>=SZ_X || y>=SZ_Y)
	return; // less than zero with signTrick
    unsigned int *ptr = (unsigned int*)screen->pixels;
    int lineoffset = y * (screen->pitch / 4);
    ptr[lineoffset + x] = color;
}


void cls() {
    int i,yofs,dy,j,ofs;
    for (i = 0, yofs = 0, dy = 0; i < SZ_Y; i++, yofs += screen->pitch / 4)
	for (j = 0, ofs = yofs; j < SZ_X; j++, ofs++, dy++)
	    ((unsigned int*)screen->pixels)[ofs] = 0;
}


// GUI vars global. code is getting messy and ugly.
// this is a hack not an application, (so) this is not an issue.
int boxSpace=(int)(400./nMass);

int scale_k=60;

void render() {   
    // Lock surface if needed
    if (SDL_LockSurface(screen) < 0) 
	return;

    // Ask SDL for the time in milliseconds
    int tick = SDL_GetTicks();

    cls(); // instade of costly CLS - we can erase blocks if needed

    if (!pause) {
	for (int n=0;n<300;++n) {
	    game_tick();
	}

	updateTract();
	tGlobal+=0.01;
    }

    SDL_Rect r;

    // draw mid-line
    r.y=SZ_Y/2+8;
    r.x=52;
    r.h=2;r.w=400;
    SDL_FillRect(screen,&r,0x00b04080);

    if ((unsigned)mouseColumn<nMass){
	r.x=50+mouseColumn*boxSpace;
	r.y=30;r.w=boxSpace-4;r.h=SZ_Y-120;
	SDL_FillRect(screen,&r,0x00ffccdd);
    }

    for (int i= 1 ; i<nMass-1; ++i ) {
	r.x=50 + i*boxSpace;
	r.y=SZ_Y/2+ (int) (tract[i].x * sc);
	r.h=16;r.w=boxSpace-4;
	if (skipMass[i])
	    SDL_FillRect(screen,&r,
		    0x00f00090+(((int)(tGlobal*100))<<0));
	else
	    SDL_FillRect(screen,&r,
		    0x0090f000+(abs((int)(tract[i].v*180))<<16));
    }

    // draw tract bounds
    r.x=450;
    r.y=SZ_Y/2-100;
    r.h=200;r.w=4;
    SDL_FillRect(screen,&r,0x00ffc000);
    r.x=46;
    SDL_FillRect(screen,&r,0x0000c0ff);

    // draw Spring k's
    for (int i= 1 ; i<nMass-1; ++i ) {
	r.x=50 + i*boxSpace;
	r.y=SZ_Y - scale_k;
	r.h=(int)(scale_k*tract[i].k); r.w=boxSpace-4;
	SDL_FillRect(screen,&r,
		    0x0000a070);
    }

    // Unlock if needed
    if (SDL_MUSTLOCK(screen)) 
	SDL_UnlockSurface(screen);

    // Tell SDL to update the whole screen
    SDL_UpdateRect(screen, 0, 0, SZ_X, SZ_Y);    

}

struct audioGen {
    virtual void audioLoop(Sint16*stream,int len)=0;
};

struct audioOsc : audioGen {
    double _phase;
    double _freq;
    audioOsc(double freq_=110.) {
	_phase=0.;
	_freq=freq_;
    }
};

struct sinOsc : audioOsc {
    sinOsc():audioOsc() {}
    sinOsc(double d):audioOsc(d) {}
    void audioLoop(Sint16 *stream,int len) {
	for (int i=0;i<len;) { 
	    stream[i++] = (int)(10000*sin(_phase));
	    _phase += _freq*PI_2*SF; // pitch
	}
    }
};

struct trainOsc : audioOsc {
    trainOsc():audioOsc() {}
    // do we have have to declare such ctor forwards?!
    trainOsc(double d):audioOsc(d) { _phase=PI_2;} // first click
    void audioLoop(Sint16 *p,int len) {
	for (int i=0;i<len;++i) { 
	    p[i]=0;
	    if (_phase >= PI_2) {
		p[i] = 12000;
		_phase -=PI_2;
	    } // pulse-train
	    _phase += _freq*PI_2*SF;
	}
    }
};

struct noiseOsc : audioOsc {
    noiseOsc():audioOsc() {}

//TODO: some cfg, testing (spec.. etc)
    void audioLoop(Sint16 *p,int len) {
	for (;--len;) { 
	    *++p=rand()&(0x8000|0x3fff);
	}
    }
};

struct naiveTriOsc : audioOsc {
    naiveTriOsc() : audioOsc() {}
    naiveTriOsc(double d):audioOsc(d) {}
    void audioLoop(Sint16 *p,int len) {
	for (int i=0;i<len;++i) { 
	    p[i]=(Sint16)(10000*_phase/PI_2);
	    if (_phase >= PI_2) {
		_phase -=PI_2;
	    } // Triangle, with aliasing
	    _phase += _freq*PI*SF; // PI or PI*2?
	}
    }
};

struct combFilter : audioGen {
    int _samplesDelay;
    Sint16 *_lastBuff;

    combFilter(int samplesDelay_) {
	_samplesDelay=samplesDelay_;
	if (_samplesDelay > aBuffSz){
	    fprintf(stderr,"combFilter!! crying!\n");
	    exit(-1);
	}

	_lastBuff=new Sint16[aBuffSz];
    // static allocation to allow dynamic delay changes
	for (int i=0;i<aBuffSz;++i)
	    _lastBuff[i]=0;
    }
    combFilter(double delaySec) {
	_samplesDelay =(int)( SR*delaySec);
	if (_samplesDelay > aBuffSz){
	    fprintf(stderr,"combFilter!! crying!\n");
	    exit(-1);
	}
	_lastBuff=new Sint16[aBuffSz];
	for (int i=0;i<aBuffSz;++i)
	    _lastBuff[i]=0;
    }
    virtual ~combFilter() {
	delete[]_lastBuff;
	_lastBuff=0;
    }

    void audioLoop(Sint16 *p,int len) {
	for (int i=0;i<_samplesDelay;++i) { 
	    p[i] += (int)(0.997*_lastBuff[i]);
	}
	for (int i=_samplesDelay;i<len;++i) { 
	    p[i] += (int)(0.997*p[i-_samplesDelay]);
	}
	for (int i=len-_samplesDelay;i<len;++i) { 
	    _lastBuff[i-len+_samplesDelay] = p[i];
	}
    }
};

// @:-  The Kelly-Lochbaum Impl:
struct KL_Filter : audioGen {
// Kelly-Lochbaum tract
// direct 4-mul lattice implementation

/*
 *
 * k = reflection coefficient
 *
 *      A(i) - A(i+1) 
 * k = --------------
 *      A(i) + A(i+1)
 *
 *  (A = Area of tube section)
 *
 *             ..-->|z^-1|----->*-------+->|z^-1|--->  f
 *                          v 1+k       ^
 *                          |           |
 * glottis ...              * k         |          .... lips
 *                          |        -k *
 *                          |           |
 *                          v     1-k   |
 *         ...----<|z^-1|<--+--<---*----|-<|z^-1|<---  b
 *
 *
 */

// waveGuide, forward and backward
#define _KL_MAX 64
    double _f[_KL_MAX]; // glottis to lips
    double _b[_KL_MAX]; // lips to glottis
    double _last[_KL_MAX]; // backup for next loop
    int _lasti;
    int _n_sections;
    double _k[_KL_MAX]; // reflection coefs

    KL_Filter() {
	this->clear();
	_n_sections=0;
    }

    void clear() {
	for (int i=0;i<_KL_MAX;++i) {
	    _f[i]=_b[i]=_last[i]=0.;
	}
	_lasti=0;

    }

// @:+ reflection from area
    void updateReflectionCoefs(double *area,int len) {
	_n_sections=len;
	if (len>_KL_MAX) {
	    fprintf(stderr,"KellyLochBaum!! too many sections given!\n");
	    exit(-1);
	}
	_k[0]=0.9997; // (closed) glottis reflection
	for (int i=1;i<len;++i) {
//	    printf("AREAs: %f , %f ",area[i-1],area[i]);
//	    printf("-, +: %f , %f\n ",area[i-1]-area[i],area[i-1]+area[i]);

	    if (0.==area[i-1] && 0.==area[i])
		_k[i]=0.;
	    else {
		double sum,diff;
		diff= area[i-1]-area[i];
		sum= area[i-1]+area[i];
		_k[i] = diff/sum;
	    }
	}

	_k[len]=-0.9997; // (open) lips reflection
    }

    void printCoefs() {
	/// dbg
	printf("KL::refl::");
	for (int i=0;i<_n_sections+1;++i) {
	    printf("(%d)%f,",i,_k[i]);
	}
	printf("\n");
	fflush(stdout);
    }


// we are excited by a 'glottis', and output a 'lips':= f+b
//
// NOTE: this is not an 'efficient' cyclic buff impl!
// but a pseudocodic memcpy impl.

// @:+ The KL ladder filter
    void audioLoop(Sint16 *p,int len) {
	int i;double out;
	int readIndex;
	i=_lasti;
	for (;len;--len,++p) {
	    _f[0]= *p; // read exciting new data
	    // should we normalise to 0..1?

	    *p = (Sint16)_last[i]; // feed prev proccessed samples
	    ++i; i %= _n_sections;
	    // TODO: should we do a double sized array trick?

	    // filter:
	    for (i=1;i<_n_sections;++i) {
		_f[i] = (1.+_k[i])*_f[i-1] -_k[i]*_b[i];
		_b[i-1] = (1.-_k[i])*_b[i] +_k[i]*_f[i-1];
	    }
	    //_b[i-1] = _k[i]*_f[i-1]; // k is 1 longer, with the lip coef
	    _b[_n_sections-1] = _k[_n_sections]*_f[_n_sections-1]; 
	    // last time for 'b' is eval'd one behind

	    // lip presure sum
	    out = _f[_n_sections-1] ;
		// + _b[_n_sections-1];
	    _last[i] = out; // to be played _n_sections samples later
	}
	_lasti=i;
    }

};


// Audio net:
audioOsc *osc0;

combFilter filt0(54);
//combFilter filt0(0.01);
KL_Filter tr0;

//void oscOrder = {
//   trainOsc,sinOsc,naiveTriOsc };
// we can't save list of ctor's in C++ can we?

audioOsc *oscs[]={new trainOsc(110.0),
    new noiseOsc,
    new sinOsc(101.0),
    new naiveTriOsc(111.0) }; // should they share freq?

int oscInd=0;

void init_audioNet() {
    osc0=oscs[oscInd];
}

void changeOsc() {
    ++oscInd;
    oscInd %= 4; /// sizeof(oscs); // doesn't work with sizeof XXX chk
    osc0=oscs[oscInd];
}


// we get a clear'd buffer for some reason.
// in a lower level system we prob could do circular buff
// => lower latency
// (expr in Win32)
void audioCall(void *unused, Uint8 *stream, int len) {
    Sint16 *p=(Sint16*)stream; len /=2; // /2 for 16 bit
    osc0->audioLoop(p,len);
//    filt0.audioLoop(p,len);
    tr0.audioLoop(p,len);
}


// Mouse UI state:
#define M_SELECT_MASS 1
int mouseState = M_SELECT_MASS;
struct Mouse {
    int x,y;
    int but0; // buttons down
};
Mouse mouse;

void init_gui() {
    boxSpace=(int)(400./nMass);
}

void init_audioVars() {
    SF = 1./SR;
    double speedOfSound = 34500. ; // cm/sec
    cmPerSample = speedOfSound*SF;
    
    tractLen=17.;
    //int tractSamplesDelay=tractLen/cmPerSample;
    int tractSamplesDelay=(int)
	(tractLen*(double)SR/speedOfSound);
    printf("Every sample delays: %f cm\n",cmPerSample);
    printf("A one second delay, would mean, SR delay == %d \n",SR);
    printf("Meaning a distance of %f cm\n",cmPerSample*SR);
    printf("So a %f cm tube, means z^-%d\n",tractLen,tractSamplesDelay);
    filt0._samplesDelay=tractSamplesDelay;

    plusSize=limit_size?minSize:0.;
}

void updateTract() {
    double area[nMass-2];
    for (int i=0;i<nMass-2;++i) {
	area[i]=plusSize+10*fabs(tract[i+1].x);
    }

    tr0.updateReflectionCoefs(area,nMass-2);
}

// Entry point
int main(int argc, char *argv[]) {
    // Initialize SDL's subsystems
    if ( SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO| SDL_INIT_JOYSTICK ) < 0 ) {
	fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
	exit(1);
    }
    
    init_masses();
    init_audioVars();
    init_audioNet();
    init_gui();

    // Register SDL_Quit to be called at exit; makes sure things are
    // cleaned up when we quit.
    atexit(SDL_Quit);
    int n=SDL_NumJoysticks();
    SDL_Joystick *joystick;

    printf("%i joysticks were found.\n", n);
    for(int i=0; i < SDL_NumJoysticks(); i++ ) {
	printf("    %s\n", SDL_JoystickName(i));

	SDL_JoystickEventState(SDL_ENABLE);
	joystick = SDL_JoystickOpen(0);      
	printf("Opening Joystick 0\n");
	break;
    }
    fflush(stdout);


    // Attempt to create a window with 32bit pixels.
    screen = SDL_SetVideoMode(SZ_X, SZ_Y, 32, SDL_HWSURFACE);

    // If we fail, return error.
    if ( screen == NULL )     {
	fprintf(stderr, "Unable to set video: %s\n", SDL_GetError());
	exit(1);
    }
    // AUDIO setup
    SDL_AudioSpec desired, obtained;
    desired.freq = SR;

    desired.format = AUDIO_S16;
    desired.samples = aBuffSz; // if too small, we miss cycles
    // ^ is this bytes? (yes! so must be even)
    desired.channels = 1;

    desired.callback=audioCall;
    desired.userdata=NULL;


//    if ( SDL_OpenAudio(&desired, NULL) < 0 ){ // force desired!
    if ( SDL_OpenAudio(&desired, &obtained) < 0 ) {
	fprintf(stderr, "AudioMixer, Unable to open audio: %s\n"
		, SDL_GetError());
	exit(1);    
    }
    if (obtained.format!=AUDIO_S16) {
	fprintf(stderr, "Didn't get Signed16bit!");
	exit(1);
    }   
    if (obtained.freq!=desired.freq) {
	fprintf(stderr, "Didn't get desired SR! [%d], what to do!?",
		desired.freq);
	exit(1);
    }   

    // -- End (AUDIO setup)

    SDL_WM_SetCaption("Mass!Tst",0);

    bool keyDown[SDLK_LAST];
    for (int i=0;i<SDLK_LAST;i++){
	keyDown[i]=false;
    }

    updateTract();
    tr0.printCoefs();

    SDL_PauseAudio(0); // run audio

    // Main loop: loop forever.
    while (1) {
	// Render stuff
	render();
	updateTract();
	for (int i=0;i<nMass;++i)
    	mousePull[i]=0;

	if (mouse.but0==1) { // down
	    if ((unsigned)mouseColumn<nMass-1 && mouseColumn ) {
		if (mouse.y < SZ_Y - scale_k ) {
		    if (skipMass[mouseColumn]) goto skipMass;
		    tract[mouseColumn].x=(mouse.y-SZ_Y/2)/sc;
		    tract[mouseColumn].v=0;
		} else { //k
		    tract[mouseColumn].k=(mouse.y-SZ_Y+60)/
			(double)scale_k;
		}
	    }
	}
skipMass:
	if (mouse.but0==2) { // down
	    if ((unsigned)mouseColumn<nMass-1 && mouseColumn) {
		mousePull[mouseColumn]=(mouse.y-SZ_Y/2)/sc;
	    }
	}

	if (keyDown[SDLK_c]) {
	    static int cCnt=11;
	    if (cCnt==0) {
		--nMass;
		if (2==nMass)
		    nMass=3;
		tract[nMass].v=tract[nMass].x=0.;
		init_gui();
		updateTract();
		//init_audioVars();
		cCnt=11;
	    }
	    cCnt--;
	}
	if (keyDown[SDLK_v]) {	
	    static int vCnt=11;
	    if (vCnt==0) {
		++nMass;
		if (65==nMass)
		    nMass=64;
		init_gui();
		updateTract();
		//init_audioVars();
		vCnt=11;
	    }
	    vCnt--;
	}
	if (keyDown[SDLK_s]){
	    if (filt0._samplesDelay<aBuffSz)
		filt0._samplesDelay++;
	}
	if (keyDown[SDLK_a]){
	    if (filt0._samplesDelay!=0)
		filt0._samplesDelay--;
	}

	if (keyDown[SDLK_f]){
	    osc0->_freq += 0.1;
	}
	if (keyDown[SDLK_d]){
	    osc0->_freq -= 0.1;
	}

//	SDL_Delay(10); // sleep, not to eat CPU
	// Poll for events, and handle the ones we care about.
	SDL_Event event;
	while (SDL_PollEvent(&event)) {
//	while (SDL_WaitEvent(&event)) {
	    switch (event.type) {
		default:
		    break;
		case SDL_KEYDOWN:
		    keyDown[event.key.keysym.sym]=true;
		    switch( event.key.keysym.sym ){
			case SDLK_z:
			    damp +=0.00001;
			    if (damp>1.) damp=1.;
			    break;
			case SDLK_x:
			    damp -=0.00001;
			    if (damp<0.) damp=0.;
			    break;
			case SDLK_c:
			    break;
			case SDLK_v:
			    break;
			case SDLK_p:
			    pause=!pause;
			    break;
			case SDLK_u:
			    updateTract();
			    tr0.printCoefs();
			    break;
			case SDLK_y:
			    tr0.clear();
			    break;
			case SDLK_e:
			    limit_size=!limit_size;
			    plusSize=limit_size?minSize:0.;
			    break;
			case SDLK_w:
			    sc*=4.;
			    break;
			case SDLK_q:
			    sc/=4.;
			    break;
			case SDLK_o:
			    changeOsc();
			    break;
		    }
		    break;

		case SDL_KEYUP: 
		    keyDown[event.key.keysym.sym]=false;
		    switch( event.key.keysym.sym ){
			case SDLK_x:
			    break;
			case SDLK_z:
			    break;
		    }

		    // If escape is pressed, return (and thus, quit)
		    if (event.key.keysym.sym == SDLK_ESCAPE)
			goto cleanup;
		    break;

		case SDL_MOUSEMOTION:
// here we go, dirty UI state machine coding.
		    mouse.x=event.motion.x;
		    mouse.y=event.motion.y;
		    switch (mouseState) {
			case M_SELECT_MASS:
			    int y=event.motion.x;
			    if (y<50) mouseColumn=-1;
			    else mouseColumn=(y-50)/boxSpace;
			    break;               
		    }
		    break;
		case SDL_MOUSEBUTTONDOWN:
		    switch (event.motion.state) {
			case SDL_BUTTON_LEFT:
			    mouse.but0=1;
			    if ((unsigned)mouseColumn < 
				    nMass-1 && mouseColumn) {
				if (mouse.y < SZ_Y - scale_k ) {
				    tract[mouseColumn].x=(mouse.y-SZ_Y/2)/sc;
				    tract[mouseColumn].v=0;
				}
			    }
			    break;
			case   SDL_BUTTON_RIGHT:
			    mouse.but0=2;
			    skipMass[mouseColumn]=!skipMass[mouseColumn];
			    break;
			    if ((unsigned)mouseColumn<nMass-1 && mouseColumn) {
				mousePull[mouseColumn]=(mouse.y-SZ_Y/2)/sc;
			    }
			    break;
		    }
		    break;
		case SDL_MOUSEBUTTONUP:
		    mouse.but0=0;
		    if ((unsigned)mouseColumn<nMass) {
			tract[mouseColumn].v=(mouse.y-event.motion.y)/sc;
		    }
		    break;

		case SDL_QUIT:
		    goto cleanup;
	    }
	}
    }
cleanup:

    return 0;
}



