#include <Bela.h>
#include <cmath>
#include <algorithm>
#include <string.h>

float gInterval = 1;
float gSecondsElapsed = 0;
int gCount = 0;

// AUDIO
#define gAUDIOFRAMES 128
int gAudioChannelNum; // number of audio channels to iterate over
float gBlock[gAUDIOFRAMES*2]; //interleaved AudioBlock to be used in render.



// TEST
float test[gAUDIOFRAMES*2];

// ------------------------- SETUP -------------------------
bool setup(BelaContext *context, void *userData)
{
    
    // Print a string
	//rt_printf("Audioframes = %d \n",context->audioFrames);
	//rt_printf("sizeof(context->audioIn)= %d \n", sizeof(context->audioIn));


	// Set number of audio channels to minimum of I/O
    gAudioChannelNum = std::min(context->audioInChannels, context->audioOutChannels);
    //set AudioFrame 
	return true;
}


// ------------------------- RENDER -------------------------
void render(BelaContext *context, void *userData)
{
	
	// tilsvarende til audioRead, men 1 samlet operation. gBlock er interleavet LR
	memcpy(gBlock, context->audioIn,sizeof(context->audioIn)*context->audioFrames*2);
	
	
	/*
	for(unsigned int n = 0; n < context->audioFrames; n++) {
	
		for(unsigned int ch = 0; ch < gAudioChannelNum; ch++){
			test[2*n+ch+1] = audioRead(context,n,ch); // Off by one at the start samp(1) == samples (0)
			audioWrite(context, n, ch, audioRead(context, n, ch));
		} //end inner audio loop
		
		// Increment a counter on every frame
		gCount++;
		
		// Print a message every second
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
	
	
	memcpy(context->audioOut,gBlock,sizeof(context->audioIn)*context->audioFrames*2);
	
}



// ------------------------- CLEANUP -------------------------
void cleanup(BelaContext *context, void *userData)
{

}


/**
\example print/render.cpp

Printing to the console
-----------------------

This example demonstrates how to print to the console. When working within the audio thread
use the function rt_printf(). This has the same functionality as printf() but is safe to call
from the audio thread. However, make sure to not make too many calls to this function within a
render loop as this may overload the CPU and/or stall communication with the board.
In the render() function above a counter is implemented in order to only print to the console
after a specified interval has passed. A counter variable is used to keep track of the amount
of samples elapsed after starting the program.
The usage of rt_printf() is identical to printf(): http://en.cppreference.com/w/cpp/io/c/fprintf

*/
