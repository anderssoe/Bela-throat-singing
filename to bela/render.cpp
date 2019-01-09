#include <Bela.h>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <biquad.h> // simple biquad filters
#include <math_neon.h>
#include "aux_functions.h"



// AUDIO VARIABLES

// Setup by user
#define gAUDIOFRAMES 128
#define FFTSIZE 1024
float gGainClean_dB =00000;
float gGainFX_dB = 0000;
int gAudioChannelInput = 0; //which input channel to use. 0 = left, 1 = right



// Internal only
int gAudioChannelNum; // number of audio channels to iterate over
float gBlock[2*FFTSIZE]; // AudioBlock to be used in render. first FFTSIZE is for input buffering, last FFTSIZE for output
float gBlockClean[gAUDIOFRAMES]; // Block to be used for clean signal
float gOutPutBuffer[gAUDIOFRAMES];
float gGainClean_lin = pow(10, gGainClean_dB / 20.0);
float gGainFX_lin = pow(10, gGainFX_dB / 20.0);
int gStatesTotal = FFTSIZE/gAUDIOFRAMES /2; // Amount of states needed
int gState = 0; // Current state
float gWindow[FFTSIZE];
float gOverlapBuffer[FFTSIZE/2] = {};


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
	// TEST

	// Set number of audio channels to minimum of I/O
    gAudioChannelNum = std::min(context->audioInChannels, context->audioOutChannels);
 	//Create window
	hann_window(gWindow, FFTSIZE,true);


    // Print strings
    rt_printf("Using %s input\n", gAudioChannelInput ? "Right":"Left");
	rt_printf("Audioframes = %d samples per input, FFTsize = %d \n",context->audioFrames,FFTSIZE);
	rt_printf("Bela buffers are %s\n", (context->flags & BELA_FLAG_INTERLEAVED) ? "interleaved":"not interleaved");
	rt_printf("Clean audio is at %f dB = %f lin \n",gGainClean_dB,gGainClean_lin);
	rt_printf("FX audio is at %f dB = %f lin \n",gGainFX_dB,gGainFX_lin);

	return true;
} //END setup


// ------------------------- RENDER -------------------------
void render(BelaContext *context, void *userData)
{
	// Copies input channel to processing block and a block holding latest clean signal
	memcpy(&gBlock[FFTSIZE/2 + gState * gAUDIOFRAMES], &context->audioIn[gAudioChannelInput * context->audioFrames],sizeof(context->audioIn)*context->audioFrames);
	memcpy(gBlockClean, &context->audioIn[gAudioChannelInput * context->audioFrames],sizeof(context->audioIn)*context->audioFrames);



	// Mix clean and processed fx signal
	for(unsigned int idx = 0; idx <gAUDIOFRAMES;idx++){

		gOutPutBuffer[idx] =gGainFX_lin*gBlock[FFTSIZE + gState*gAUDIOFRAMES + idx] + gGainClean_lin* gBlockClean[idx];
	}

	// Send total signal to outputs
	memcpy(context->audioOut,gOutPutBuffer,sizeof(context->audioIn)*context->audioFrames);	//Copy output block to left output channel
	memcpy(&context->audioOut[1 * context->audioFrames],gOutPutBuffer,sizeof(context->audioIn)*context->audioFrames);  //Copy output block to right output channel


	// Increment state
	gState++;

	if(gState>=gStatesTotal){
		gState = 0;

		// Copy full block to processing part
		memcpy(&gBlock[FFTSIZE], &gBlock[0],sizeof(context->audioIn)*FFTSIZE);

		//Move newest 4 blocks, to be ready for next iteration
		memcpy(gBlock,&gBlock[FFTSIZE/2],sizeof(context->audioIn)*FFTSIZE/2);

		// Analysis window
		for(unsigned int idx = 0;idx<FFTSIZE;idx++){
			gBlock[FFTSIZE + idx] = gBlock[FFTSIZE + idx]*gWindow[idx];
		}


		// PROCESS HERE


		// Synthesis window
		for(unsigned int idx = 0;idx<FFTSIZE;idx++){
			gBlock[FFTSIZE+idx] *=gWindow[idx];
		}

		//overlap-add
		for(int idx = 0;idx<FFTSIZE/2;idx++){
			gBlock[FFTSIZE + idx] +=gOverlapBuffer[idx];
		}
		// save new overlap
		memcpy(gOverlapBuffer,&gBlock[FFTSIZE + FFTSIZE/2 ],sizeof(context->audioIn)*FFTSIZE/2);

	}
} // END render

// ------------------------- CLEANUP -------------------------
void cleanup(BelaContext *context, void *userData)
{

}
