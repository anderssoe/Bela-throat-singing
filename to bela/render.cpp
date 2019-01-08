#include <Bela.h>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <biquad.h> // simple biquad filters
#include <math_neon.h>

// AUDIO VARIABLES

// Setup by user
#define gAUDIOFRAMES 1024
#define FFTSIZE 128
float gGainClean_dB = 0;
float gGainFX_dB = 0;
int gAudioChannelInput = 0; //which input channel to use. 0 = left, 1 = right



// Internal only
int gAudioChannelNum; // number of audio channels to iterate over
float gBlock[2*FFTSIZE]; // AudioBlock to be used in render. first FFTSIZE is for input buffering, last FFTSIZE for output
float gBlockClean[gAUDIOFRAMES]; // Block to be used for clean signal
float gOutPutBuffer[gAUDIOFRAMES];
float gGainClean_lin = pow(10, gGainClean_dB / 20.0);
float gGainFX_lin = pow(10, gGainFX_dB / 20.0);
int gStatesTotal = FFTSIZE/gAUDIOFRAMES; // Amount of states needed 
int gState = 0; // Current state



// TEST VARIABLES

// Variables used to print every gInterval during render
float gInterval = 1; 
float gSecondsElapsed = 0;
int gCount = 0;


// ------------------------- USER SETTINGS -------------------------
void Bela_userSettings(BelaInitSettings *settings) 
{
	settings->interleave = 0; //Sets if the audio buffers are interleaved or not. 0 = not interleaved.
	settings->analogOutputsPersist = 0; //makes sure analog outs do not persist across frames. Reduces memory usage
}


// ------------------------- SETUP -------------------------
bool setup(BelaContext *context, void *userData)
{
    // Print strings
    rt_printf("Using %s input\n", gAudioChannelInput ? "Right":"Left");
	rt_printf("Audioframes = %d samples per input \n",context->audioFrames);
	rt_printf("Bela buffers are %s\n", (context->flags & BELA_FLAG_INTERLEAVED) ? "interleaved":"not interleaved");
	rt_printf("Clean audio is at %f dB = %f lin \n",gGainClean_dB,gGainClean_lin);
	rt_printf("FX audio is at %f dB = %f lin \n",gGainFX_dB,gGainFX_lin);
	// Set number of audio channels to minimum of I/O
    gAudioChannelNum = std::min(context->audioInChannels, context->audioOutChannels);
 
	return true;
} //END setup


// ------------------------- RENDER -------------------------
void render(BelaContext *context, void *userData)
{
	// Copies input channel to processing block and a block holding latest clean signal
	memcpy(&gBlock[gState * gAUDIOFRAMES], &context->audioIn[gAudioChannelInput * context->audioFrames],sizeof(context->audioIn)*context->audioFrames);
	memcpy(gBlockClean, &context->audioIn[gAudioChannelInput * context->audioFrames],sizeof(context->audioIn)*context->audioFrames);
		
	
	
	// Mix clean and processed fx signal
	for(unsigned int idx = 0; idx <gAUDIOFRAMES;idx++){ 
		gOutPutBuffer[idx] =gGainFX_lin*gBlock[FFTSIZE + gState*gAUDIOFRAMES + idx] + gGainClean_lin* gBlockClean[idx];
	}
	
	// Send total signal to outputs
	memcpy(context->audioOut,gOutPutBuffer,sizeof(context->audioIn)*context->audioFrames);	//Copy output block to left output channel
	memcpy(&context->audioOut[1 * context->audioFrames],gOutPutBuffer,sizeof(context->audioIn)*context->audioFrames);  //Copy output block to right output channel
	
	
	// Increment state and copy a full input buffer to processing
	gState++;
	if(gState >= gStatesTotal){
		gState = 0;
		memcpy(&gBlock[FFTSIZE], &gBlock[0],sizeof(context->audioIn)*FFTSIZE);
	}
} // END render

// ------------------------- CLEANUP -------------------------
void cleanup(BelaContext *context, void *userData)
{

}

