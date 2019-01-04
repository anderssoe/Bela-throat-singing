#include <Bela.h>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <biquad.h> // simple biquad filters
#include <math_neon.h>

// AUDIO VARIABLES
#define gAUDIOFRAMES 128
int gAudioChannelNum; // number of audio channels to iterate over
int gAudioChannelInput = 0; //which input channel to use. 0 = left, 1 = right

float gBlock[gAUDIOFRAMES]; // AudioBlock to be used in render.



// TEST
float test[gAUDIOFRAMES];

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
    rt_printf("Using input %s \n", gAudioChannelInput ? "Right":"Left");
	rt_printf("Audioframes = %d samples per input \n",context->audioFrames);
	rt_printf("Bela buffers are %s\n", (context->flags & BELA_FLAG_INTERLEAVED) ? "interleaved":"not interleaved");


	// Set number of audio channels to minimum of I/O
    gAudioChannelNum = std::min(context->audioInChannels, context->audioOutChannels);



	return true;
}


// ------------------------- RENDER -------------------------
void render(BelaContext *context, void *userData)
{

	// Copies input channel to processing block
	memcpy(gBlock, &context->audioIn[gAudioChannelInput * context->audioFrames],sizeof(context->audioIn)*context->audioFrames);





/*	for(unsigned int n = 0; n < context->audioFrames; n++) { //main audio loop if not using block

		for(unsigned int ch = 0; ch < gAudioChannelNum; ch++){



			audioWriteNI(context, n, ch, audioReadNI(context, n, ch));


		} //end inner audio loop

		// Increment a counter on every frame
		gCount++;

		// Print a message every interval
		if(gCount % (int)(context->audioSampleRate*gInterval) == 0) {
		    gSecondsElapsed += gInterval;



	    	rt_printf("memcpy: ");
		    for(int j = 0;j<context->audioFrames;j++){
				rt_printf("%f, ",gBlock[j]);
		    }
		    rt_printf("\n audioRead: ");
		    for(int j = 0;j<context->audioFrames;j++){
				rt_printf("%f, ",test[j]);
		    }
		}



	} //end outer audio loop
*/

	memcpy(context->audioOut,gBlock,sizeof(context->audioIn)*context->audioFrames);	//Copy block to left output channel
	memcpy(&context->audioOut[1 * context->audioFrames],gBlock,sizeof(context->audioIn)*context->audioFrames);  //Copy block to right output channel

}



// ------------------------- CLEANUP -------------------------
void cleanup(BelaContext *context, void *userData)
{

}
