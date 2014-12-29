This is a GUI that manages DSP analysis functions for wav-files (especially speech signals). Two functions (plotps.m and spect.m) are included for starters. You may write your own functions and integrate that into the GUI without much hassle (see instructions in the accompanying Readme.txt file). Additional features like the snipper lets you trim the time series and save it as a separate wav-file. The GUI is a great tool for instructors in a DSP course and DSP researchers alike!

If you like the GUI, please comment below; otherwise (e.g., if you find a suspicious bug), send an email to me (weechiat@gmail.com). :)

Essential files: spsnip_gui.p, spsnip_guicb.p, spsnip_config.txt
Accompanying files: plotps.m, spect.m, spectmesh.m, spectgray.m, pitchwatch.m, Hum.wav, Readme.txt

** Updates **
[25/07/06] Added spectmesh.m (mesh-like spectrogram) and spectgray.m (grayscale spectrogram).
[03/05/07] Added pitchwatch.m.
[18/11/07] Drag-drop cursors replacing slider-bars.

** General instructions **
1. To launch, type "spsnip_gui" in the workspace.
2. To use the snipper, trim the time-series with the drag-drop red cursors to the desired start- and end-points. Then, click on the "Save..." button to save the portion of the time-series as a wav-file.

** Additional features **
The GUI will register its auxiliary functions (listed in spsnip_config.txt) every time it is launched. These functions have the function syntax:
> function_name(x,Ts)
where x is the input (vector) time-series and Ts (in sec) is the sampling time. If you intend for your function to plot a figure, do create a new (blank) figure before plotting the figure. Otherwise, you run the risk of overwriting the figure in the GUI. Finally, do post whatever interesting functions you create here at Mathworks for all to share.

** Limitations **
1. The GUI may be buggy if the loaded wav-file is stereo/multichannel.
