/*	TRACK RUPTURES 
  	***********************
  	
  	Author and copyright: 		Winnok H. De Vos
  	Date Created: 				April 06, 2013
  	Date Last Modified:			Nov 25, 2022
  	
  	This macro is conceived to track nuclear labeled cells and detect nuclear ruptures 	
  	by measuring signal intensity changes of a GFP-NLS signal. 
  	Works with single channel images as well as 
  	multichannel images with one permanent nuclear marker such as a histone fusion protein.
  	Change log:
  	+ Added multiple track correction tools
  	+ Measures two instead of selected channel
  	+ Added manual ROI addition tool
  	+ Corrected an issue with rois not showing on hyperstack 
  	+ added stardist for improved nuclear segmentation
  	+ revised ROI issues and changed tracking jargon (Frames instead of Slices)
	+ corrected the ROI assignemnt after Stardist detection to allow separation of overlapping nuclei (30/03/21)
	+ nuclear shape
	
  	Please cite the relevant paper if useful for your research:
  	J. Robijns, F. Molenberghs, T. Sieprath, T. Corne, M. Verschuuren and W. De Vos (2016). 
  	In silico synchronization reveals regulators of nuclear ruptures in lamin A/C deficient model cells. 
  	Scientific reports 6, 30325, p. 1-11. 
  	  	
  	**********************************************
*/


/*
  *********************************************
  	
  		Variable declaration 
  				
  **********************************************
*/

//	Arrays
var extent 					= newArray("This Time Point Only","All Time Points","
							+"All Previous Time Points","All Following Time Points");	
																	//	Extent of track editing			
var enhancers				= newArray("Gauss","Median","Laplace");	//	Blob enhancement algorithm
var thresholds 				= getList("threshold.methods");			//	Autothreshold methods
																	
//	Integers
var channels				= 0;									//	Number of channels
var enhance					= 50;									//	CLAHE blocksize nuclei_radiusius
var frames					= 0;									//	Number of frames
var gap						= 5;									//	Time gap used for tracking with interruption
var image_height			= 1000;									//	Image image_height
var image_width				= 1000;									//	Image image_width
var maxdisp 				= 50;									//	Maximum displacement (calibrated units)
var min_nuclei_area			= 500;									// 	Min. size of nuclei in pixels	
var nr 						= 0;									//	Number of results
var nuclei_overlap 			= 0.20;									//	Amount of overlap between nuclei tolerated for Stardist nuclei detection
var nuclei_probability		= 0.35;									//	Minimal probability for Stardist nuclei detection
var nuclei_radius			= 12;									//	Spatial smoothing nuclei_radiusius
var segment_channel			= 1;									//	Segmentation Channel
var smudge					= 0;									//	Temporal smudging nuclei_radiusius (max fill)
var separate				= 0.9;									//	Separate touching cells using conditional watershed (0>v>1),no watershed (v==1), simple watershed (v==0).

//	Strings 
var autothresh 				= "Triangle";							//	Autothreshold method 
var dir 					= "";									//	Directory for saving all	
var ext 					= "This Time Point Only"; 				//	Extent of track editing
var id 						= "";									//	Image ID
var nuclei_enhancer			= "Laplace";							//	Nuclei enhancement algroithm
var prefix 					= "XY00";								//	Image prefix

//	Booleans
var ai 						= false;								//	use stardist for segmentation
var norm					= true;									//	Normalize image intensity across time

/*	
 **********************************************
  	
		Macros
  		
 **********************************************
*/

macro "Setup [s]"
{
	setup();
}

macro "Setup Action Tool - C888 T5f16s"
{
	setup();
}

macro "Nuclei Detection Action Tool - C333 o3133 of033 odd33 o0a33"
{
	//	Segment objects and return putative nuclei as ROIset
	start = getTime;
	erase();
	printSettings();
	id = getImageID; 
	setBatchMode(true);
	print("Detecting Nuclei...");
	rmc = segment(id);
	setBatchMode("exit and display");
	update_rois(id);	
	print("Done:",rmc,"objects detected");	
	toggleOverlay();
	stop = getTime;
	print("Processing Time:",(stop-start)/1000,"sec");
}

macro "Track Action Tool - C999 o0e11 o2c22 C777 o5933 C333 oa644 o0144 C777 o5033 o9022 C999 ob011"
{
	//	Tracks nuclei through time based on a ROIset
	if(roiManager("count")==0)exit("No ROIset, run detection procedure first");
	start = getTime;
	print("Tracking Nuclei...");
	setBatchMode(true);
	id = getImageID; 
	nr = measure(id, segment_channel);
	tracks = trackRegions(nr);
	coded_id = colorCode(id, tracks);
	Array.getStatistics(tracks, min, trackNr);
	tracks = sample(id);
	print("Done:",trackNr,"tracks detected");
	setBatchMode("exit and display");
	stop = getTime;
	print((stop-start)/1000, "sec");
	run("Tile");
}

macro "Select Tracks Action Tool - C333 o0e11 o2c22 o5933 oa644"
{
	setBatchMode(true);
	nr 			= nResults;
	tracks 		= newArray(nr);
	labels 		= newArray(nr);
	defaults 	= newArray(nr);
	positive 	= newArray(nr);
	indices 	= newArray(0);
	if(nr==0)exit("No results found, please (re-)run analysis");
	if(isOpen("Coded")){selectWindow("Coded");close;}
	for(i=0;i<nr;i++)
	{
		tracks[i] = getResult("Track",i);
	}
	Array.getStatistics(tracks, min, max);
	for(i=0;i<max;i++)
	{
		labels[i] = d2s(i+1,0);
		defaults[i] = 1;
	}
	rows = round(sqrt(max));
	if(pow(rows,2)>=max)columns = rows;
	else columns = rows+1;
	Dialog.create("List Selected Tracks");
	Dialog.addMessage("Check tracks to list");
	Dialog.addCheckboxGroup(rows, columns, labels, defaults)
	Dialog.show;
	for(i=0;i<max;i++)
	{
		positive[i] = Dialog.getCheckbox;
	}
	for(i=nr-1;i>=0;i--)
	{
		track = getResult("Track",i);
		if(positive[track-1]==0)
		{
			indices = Array.concat(indices,i);
			IJ.deleteRows(i,i);		
		}
	}
	for(i=0;i<indices.length;i++)
	{
		roiManager("select",indices[i]);
		roiManager("delete");
	}
	run("Select None");
	id = getImageID;
	tracks = sample(id);
	coded_id = colorCode(id, tracks);
	Array.getStatistics(positive, min, max, mean);
	total_left = mean * positive.length;
	setBatchMode("exit and display");
	print("Selected" + total_left + "tracks");
}

macro "Track Editing Tool - C333 o0e11 o2c22 Cc00 o5933 C333 oa644" 
{
	// shift+click = delete(17);  click = switch (16);
	setBatchMode(true);
	nr = nResults;
	if(nr==0)exit("No results found, please (re-)run analysis");
	getCursorLoc(x, y, z, flag); 
	//toScaled(x, y);
	if(Stack.isHyperstack)
	{
		Stack.getDimensions(image_width, image_height, channels, slices, frames);
		frame = floor(z/channels)+1;
	}
	else frame = z+1;
	nn 		= nearestNeighbour(x,y,frame,1000);
	index 	= nn[0];  							// index of the selected track in the selected time point
	track 	= getResult("Track",index); 
	time 	= getResult("Frame",index);	
	print("X:",x,"- Y:",y,"- Frame:",frame,"- Track:",track);	
	if(flag == 17 || flag == 1) 				// left button = 16; shift=1; delete track
    { 
		Dialog.create("Delete Track");
		Dialog.addNumber("Delete Track",track);
		Dialog.addChoice("Delete Track in",extent,ext); 
		Dialog.show;
		track	= Dialog.getNumber;
		ext 	= Dialog.getChoice;
		if(ext == "This Time Point Only")
		{
			IJ.deleteRows(index,index);
			roiManager("select",index);
			roiManager("delete");
			roiManager("deselect");
			run("Select None");	
			print("Deleted Track",track);
   		}	
   		else
   		{
   			for(i=nr-1;i>=0;i--)
			{
				if(getResult("Track",i)==track)
				{
					if((ext == "All Previous Time Points") && (getResult("Frame",i)<=time))	
					{IJ.deleteRows(i,i);roiManager("select",i);roiManager("delete");}
					else if ((ext == "All Following Time Points") && (getResult("Frame",i)>=time))
					{IJ.deleteRows(i,i);roiManager("select",i);roiManager("delete");}		
					else if (ext == "All Time Points")	
					{IJ.deleteRows(i,i);roiManager("select",i);roiManager("delete");}	
   				}
			}
			roiManager("deselect");
			run("Select None");	
			print("Deleted Track",track,"in",ext);
   		}	
		toggleOverlay();
		toggleOverlay();
    } 
	else if(flag == 16 || flag == 8)	// switch track
    {
		Dialog.create("Change Track");
		Dialog.addNumber("Change Track "+track+" to:",0);
		Dialog.addChoice("Change Track in",extent,ext); 
		Dialog.show;
		newtrack = Dialog.getNumber;
		ext = Dialog.getChoice;
		if(ext == "This Time Point Only")
		{
			setResult("Track",index,newtrack);
			updateResults();
			print("Changed Track",track,"to",newtrack);
   		}	
		else 
		{
			for(i=nr-1;i>=0;i--)
			{
				if(getResult("Track",i)==track)
				{
					if((ext == "All Previous Time Points") && (getResult("Frame",i)<=time))setResult("Track",i,newtrack);		
					else if ((ext == "All Following Time Points") && (getResult("Frame",i)>=time))setResult("Track",i,newtrack);		
					else if (ext == "All Time Points")setResult("Track",i,newtrack);		
				}		
			}
			updateResults();
			print("Changed Complete Track",track,"to",newtrack,"in",ext);
		}
		toggleOverlay();
		toggleOverlay();
    } 
    setBatchMode(false);
}

macro "Use Prev/Next Tool - C333 o0e11 o2c22 o5933 oa644 Cc00 oc144" 
{
	//setBatchMode(true);
	// shift+click = prev(17); click = next (16);
	nr = nResults; 
	if(nr==0)exit("No results found, please (re-)run analysis");
	im_id = getImageID;
	getCursorLoc(x, y, z, flag);
	//toScaled(x, y);
	if(Stack.isHyperstack)
	{
		Stack.getDimensions(image_width, image_height, im_channels, im_slices, im_frames);
		frame = floor(z/im_channels)+1;
	}
	else frame = z+1;
	print("X:",x,"- Y:",y,"- Frame:",frame);	
	if(flag == 17 || flag == 1) // prev roi
    { 
		prevframe = frame-1;
		nn = nearestNeighbour(x,y,prevframe,1000);
		index = nn[0]; 
		track = getResult("Track",index);
		selectImage(im_id);
		roiManager("select",index);
		Roi.getPosition(c, s, f);
		Stack.setFrame(frame);
		roiManager("Add");
		roiManager("select",roiManager("count")-1);
		Roi.setPosition(c, s, f+1);
		roiManager("Remove Channel Info");
		setResult("X",nr,x); 
		setResult("Y",nr,y); 
		setResult("Track",nr,track); 
		setResult("Frame",nr,frame); 
		updateResults();
		roiManager("deselect");
		run("Select None");
		toggleOverlay();
		toggleOverlay();
    }
   	else if(flag == 16 || flag == 8)	// next roi
    {
		nextframe = frame+1;
		nn = nearestNeighbour(x,y,nextframe,1000);
		index = nn[0]; 
		track = getResult("Track",index); 
		selectImage(id);
		roiManager("select",index);
		Roi.getPosition(c, s, f);
		Stack.setFrame(frame);
		roiManager("Add");
		roiManager("select",roiManager("count")-1);
		Roi.setPosition(c, s, f-1);
		setResult("X",nr,x); 
		setResult("Y",nr,y); 
		setResult("Track",nr,track); 
		setResult("Frame",nr,frame); 
		updateResults();
		selectImage(id);
		roiManager("deselect");
		run("Select None");
		toggleOverlay();
		toggleOverlay();
    } 
    setBatchMode(false);
}

macro "Manually Add Track ROI Action Tool - C888 T5f16m"
{
	add();
}

macro "Manually Add Track ROI [m]"
{
	add();
}

macro "Recalculate Action Tool - C888 T5f16r"
{
	setBatchMode(true);
	if(isOpen("Coded")){selectWindow("Coded");close;}
	if(isOpen("TrackData")){selectWindow("TrackData");run("Close");}
	id = getImageID;
	nr = nResults; 
	rmc = roiManager("count"); 
	tracks = sample(id);
	listTracks();
	coded_id = colorCode(id, tracks);
	setBatchMode("exit and display");
}

macro "Toggle Overlay [t]"
{
	toggleOverlay();
}

macro "Toggle Overlay Action Tool - C888 T5f16t"
{
	toggleOverlay();
}

macro "Save and Close All Action Tool - C888 T5f16c"
{
	dir = getDirectory("");
	prefix = getString("Prefix", prefix);
	selectWindow("Results");
	saveAs("Results", dir+prefix+"_Results.xls");
	run("Close");
	selectWindow("TrackData");
	saveAs("Results", dir+prefix+"_TrackData.txt");
	run("Close");
	if(isOpen("Coded"))
	{
		selectWindow("Coded");
		saveAs("Tiff", dir+prefix+"_Coded.tif");
		close();
	}
	selectWindow("Log");
	saveAs("Text", dir+prefix+"_Log.txt");
	run("Close");
	roiManager("deselect");
	roiManager("Save", dir+prefix+"_RoiSet.zip");
	roiManager("reset");
	erase();
}

macro "Save All Action Tool - C888 T0f16c T7f16-"
{
	dir = getDirectory("");
	prefix = getString("Prefix", prefix);
	selectWindow("Results");
	saveAs("Results", dir+prefix+"_Results.xls");
	selectWindow("TrackData");
	saveAs("Results", dir+prefix+"_TrackData.txt");
	if(isOpen("Coded"))
	{
		selectWindow("Coded");
		saveAs("Tiff", dir+prefix+"_Coded.tif");
	}
	selectWindow("Log");
	saveAs("Text", dir+prefix+"_Log.txt");
	roiManager("deselect");
	roiManager("Save", dir+prefix+"_RoiSet.zip");
}

macro "Save and Close All Action Tool - C888 T0f16c T7f16+"
{
	dir = getDirectory("");
	prefix = getString("Prefix", prefix);
	selectWindow("Results");
	saveAs("Results", dir+prefix+"_Results.xls");
	run("Close");
	selectWindow("TrackData");
	saveAs("Results", dir+prefix+"_TrackData.txt");
	run("Close");
	if(isOpen("Coded"))
	{
		selectWindow("Coded");
		saveAs("Tiff", dir+prefix+"_Coded.tif");
		close();
	}
	selectWindow("Log");
	saveAs("Text", dir+prefix+"_Log.txt");
	run("Close");
	roiManager("deselect");
	roiManager("Save", dir+prefix+"_RoiSet.zip");
	roiManager("reset");
	erase();
}

/*
	**********************************************

		Functions
		
	**********************************************
*/

function erase()
{
	run("Clear Results");
	roiManager("reset");
	print("\\Clear");
	run("Collect Garbage");
}

function setup()
{
	Dialog.create("Track Ruptures");
	Dialog.addNumber("Segmentation Channel", segment_channel);
	Dialog.setInsets(0,0,0);
	Dialog.addCheckbox("Normalize Intensity Over Time",norm);	
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("---------------------------------------------------")
	Dialog.setInsets(0,0,0);
	Dialog.addCheckbox("Use Stardist (overrules other segmentation settings)",ai);	
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Tolerated Overlap", 0, 1, nuclei_overlap);
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Detection Probability", 0, 1, nuclei_probability);							
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("---------------------------------------------------")
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Nuclei enhancement algorithm",enhancers, nuclei_enhancer);	
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Enhance local contrast", 0, 200, enhance);
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Nuclear nuclei radius", 0, 20, nuclei_radius);
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Temporal Smudge", 0, 200, smudge);
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Nuclear Min. Area", 0, 1000, min_nuclei_area);
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Autothreshold Algorithm",thresholds,autothresh);	
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Separation Threshold", 0, 1, separate);
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("---------------------------------------------------")
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Bridge gaps", 0, 200, gap);
	Dialog.setInsets(0,0,0);
	Dialog.addSlider("Maximum Displacement", 0, 10000, maxdisp);
	Dialog.show;
	segment_channel		= Dialog.getNumber(); 	
	norm				= Dialog.getCheckbox();
	ai					= Dialog.getCheckbox();	
	nuclei_overlap		= Dialog.getNumber(); 	
	nuclei_probability	= Dialog.getNumber(); 
	nuclei_enhancer 	= Dialog.getChoice();
	enhance				= Dialog.getNumber(); 
	nuclei_radius	 	= Dialog.getNumber(); 
	smudge	 			= Dialog.getNumber(); 		
	min_nuclei_area		= Dialog.getNumber(); 	
	autothresh 			= Dialog.getChoice();	
	separate			= Dialog.getNumber();	
	gap	 				= Dialog.getNumber(); 	
	maxdisp 			= Dialog.getNumber(); 	
	printSettings();
}

function printSettings()
{
	print("*******************************************");
	print("               Settings");
	print("*******************************************");
	print("- Segmentation channel =",segment_channel);
	print("- Normalize intensity =",norm);
	print("- Use Stardist =",ai);
	print("- Tolerated overlap =",nuclei_overlap);
	print("- Detection probability =",nuclei_probability);
	print("- Nuclei enhancement algorithm =",nuclei_enhancer);
	print("- CLAHE block size =",enhance);
	print("- Laplacian nuclei radius =",nuclei_radius);
	print("- Temporal smudge distance =",smudge);
	print("- Min. nuclear area =",min_nuclei_area);
	print("- Autothreshold algorithm =",autothresh);
	print("- Split touching nuclei =",separate);
	print("- Time gap distance =",gap);
	print("- Max displacement =",maxdisp);
	print("*******************************************");
}

function select(id,ch)
{
	selectImage(id); 
	run("Properties...", " unit=pixel pixel_width=1 pixel_height=1");
	Stack.getDimensions(image_width, image_height, channels, slices, frames); 
	if(frames==1)
	{
		Stack.setDimensions(channels, frames, slices);
		print("switched slices and frames >>>", slices, "frames"); 
	}
	selectImage(id); 
	if(Stack.isHyperstack && ch>0)		//	Multiple channels
	{
		run("Duplicate...", "title=Temp duplicate channels="+ch);
	} 
	else 								//	Single channel			
	{
		run("Duplicate...", "title=Temp duplicate");
	}
	temp_id = getImageID; 
	return temp_id;
}

function normalize(temp_id)
{
	// calculate global mean
	selectImage(temp_id); 
	global_mean = 0;
	for(i=1;i<=frames;i++)
	{
		setSlice(i);
		getRawStatistics(np,mean);
		global_mean+=mean;		
	}
	global_mean = global_mean/frames;
	
	// equalize mean intensity across time
	variance = 0;
	for(i=1;i<=frames;i++)
	{
		setSlice(i);
		getRawStatistics(np,mean);
		v = global_mean/mean;
		variance += pow(mean - global_mean,2);
		run("Multiply...", "value="+v+" slice");
	}
	variance = variance/frames;
	print("Global mean intensity across stack:", global_mean);
	print("Variance of intensities across stack:", variance);
	print("Intensity covariance across stack:", variance/global_mean);
	return temp_id;
}

function segment(id)
{
	setBatchMode(true);
	run("Collect Garbage");
	// select segmentation channel
	temp_id = select(id,segment_channel); 
	// normalize intensities over time
	if(norm)
	{
		temp_id = normalize(temp_id);
	}
	selectImage(temp_id); // generic title = Temp
	// Stardist
	if(ai)
	{	
		run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], "
		+"args=['input':'Temp', 'modelChoice':'Versatile (fluorescent nuclei)',"
		+"'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'100',"
		+"'probThresh':'"+nuclei_probability+"', 'nmsThresh':'"+nuclei_overlap+"', 'outputType':'ROI Manager', 'nTiles':'1', "
		+"'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'true', "
		+"'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");		
		newImage("Bin", "8-bit black", image_width, image_height, frames);
		binary_id = getImageID;  
		selectImage(binary_id);
		roiManager("Deselect");
		nr_rois = roiManager("count");
		for(r = 0; r < nr_rois; r++)
		{
			roiManager("select",r);
			run("Enlarge...", "enlarge=1");
			run("Clear","slice");
			run("Enlarge...", "enlarge=-1");
			run("Fill","slice");
		}
		roiManager("Deselect");
		run("Select None");
		roiManager("reset");
		setThreshold(1,255);
		run("Convert to Mask", "method=Default background=Dark black");
	}
	// Regular segmentation with some pre-processing
	else 
	{
		//	Enhance local contrast
		if(enhance>0)					
		{
			selectImage(temp_id);
			for(i=1;i<=frames;i++)
			{
				setSlice(i);
				run("Enhance Local Contrast (CLAHE)", "blocksize="+enhance+" histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
			}
		}		
		//	Smudge ruptures over time (GFP-NLS channel)
		if(smudge>0)						
		{
			selectImage(temp_id);
			run("Maximum 3D...", "x=0 y=0 z="+smudge);
		}	
		//	Blob enhancement and detection
		selectImage(temp_id);
		if(nuclei_enhancer == "Laplace")
		{
			run("FeatureJ Laplacian", "compute smoothing="+nuclei_radius);
			selectImage("Temp Laplacian");
			binary_id = getImageID;
			selectImage(binary_id);
			setAutoThreshold(autothresh+" stack");
			setOption("BlackBackground", true);
			run("Convert to Mask", "method="+autothresh+" background=Light black");
			run("Fill Holes","stack");
		}
		else
		{
			run("Duplicate...","title = Bin duplicate");
			binary_id = getImageID;
			selectImage(binary_id);
			if(nuclei_enhancer == "Gauss")run("Gaussian Blur...", "sigma="+nuclei_radius+" stack");
			else if(nuclei_enhancer == "Median")run("Median...", "radius="+nuclei_radius+" stack");
			setAutoThreshold(autothresh+" dark stack");
			setOption("BlackBackground", true);
			run("Convert to Mask", "method="+autothresh+" background=Dark black");
			run("Fill Holes","stack");
		}		
		//	Separate touching nuclei
		if(separate==0)run("Watershed","stack");
		else if (separate<1)binary_id = conditionalWatershed(binary_id,temp_id);
	}	
	//	Filter objects and clean up
	roiManager("reset");
	selectImage(binary_id);
	run("Analyze Particles...", "size="+min_nuclei_area+"-"+image_width*image_height/4+" pixel show=Nothing add exclude stack");
	rmc = roiManager("count");
	selectImage(binary_id); close;
	selectImage(temp_id); close;
	run("Collect Garbage");
	return rmc;
}

function conditionalWatershed(bin,ref)
{
	run("Colors...", "foreground=black background=black selection=yellow");
	run("Clear Results");
	roiManager("reset");
	selectImage(ref);
	title = getTitle;	
	//	Simple watershed
	selectImage(bin);
	run("Duplicate...", "title=ws duplicate");
	wid = getImageID;
	selectImage(wid); 
	run("Watershed", "stack");
	//	Detect and measure separation lines
	imageCalculator("XOR create stack",bin,wid);
	divid = getImageID;
	run("Set Measurements...", "area fit centroid median min stack redirect=["+title+"] decimal=4");
	setBatchMode("show");
	run("Analyze Particles...", "pixel display add stack");
	setBatchMode("hide");	
	//	Measure corrsponding objects and select separators with sufficient intensity drop (ratio<separate) or boundary length (<35% minor axis of fitted ellipse)
	nr = nResults;
	indices = newArray(0);
	for(i=0;i<nr;i++)
	{
		x = getResult("X",i);
		y = getResult("Y",i);
		z = getResult("Slice",i);
		selectImage(bin);
		setSlice(z);
		toUnscaled(x,y);
		doWand(x,y);
		run("Measure");
		run("Select None");
		l1 = getResult("Area",i);
		l2 = getResult("Minor",i+nr);
		m1 = getResult("Median",i);
		m2 = getResult("Median",i+nr);
		mratio = m1/m2;
		lratio = l1/l2;
		if(mratio < separate && lratio<0.35){indices = Array.concat(indices,i);}
	}
	//	Apply conditional watershed and clean up
	selectImage(bin);
	roiManager("select",indices);
	roiManager("Fill");
	selectImage(wid); close;
	selectImage(divid); close;
	run("Clear Results");
	roiManager("reset");
	run("Colors...", "foreground=white background=black selection=yellow");
	return(bin);
}

function measure(im_id, segment_channel)
{
	selectImage(im_id);
	Stack.setChannel(segment_channel);
	run("Set Measurements...", "area mean shape centroid stack redirect=None decimal=4");
	run("Clear Results");
	roiManager("deselect");
	roiManager("Measure");
	return nResults;
}

function trackRegions(nr)
{
	tracks 		= newArray(nr); 
	distances	= newArray(nr); 
	tracknr		= 0; 
	for(j = 0; j < nr; j++)
	{
		x1 			= getResult("X",j); 
		y1 			= getResult("Y",j); 
		t1 			= getResult("Frame",j);	
		a1 			= getResult("Area",j);
		if(tracks[j] == 0)
		{
			tracknr++;
			tracks[j] = tracknr;
			print(tracknr);
		} 
		if(t1<frames)
		{
			t2 = t1+1;
			nn = nearestNeighbour(x1,y1,t2,maxdisp);			//	find nearest neighbour in next time point
			while(nn[0]<0 && t2<=t1+gap) 						//	no neighbour found in next time point
			{
				t2++;
				nn = nearestNeighbour(x1,y1,t2,maxdisp); 		//	look in time span specified by time gap
			}
			index = nn[0];
			distance = nn[1];
			if(index>=0)			
			{
				if(tracks[index]==0)							// new track
				{
					tracks[index] = tracks[j];
					distances[index] = distance;
				}
				else if(distances[index]>distance)				// track exists already, choose the one with smallest distance
				{	
					tracks[index] = tracks[j];
					distances[index] = distance;
				}
			}
		}
	}
	for(j = 0;j < nr; j++)
	{
		setResult("Track",j,tracks[j]);
		setResult("Displacement",j,distances[j]);
	}
	updateResults();
	//reconnect(nr);
	return tracks;
}

function reconnect(nr)
{
	setBatchMode("Show");
	for(j = nr-1; j >= 0 ; j++)
	{
		a1 = getResult("Area",j);
		m1 = getResult("Mean",j);
		x1 = getResult("X",j);
		y1 = getResult("Y",j);
		s1 = getResult("Frame",j);
		t1 = getResult("Track",j);
		for(k = j-1; k >= 0; k++)
		{
			a2 = getResult("Area",k);
			m2 = getResult("Mean",k);
			x2 = getResult("X",k);
			y2 = getResult("Y",k);
			s2 = getResult("Frame",k);
			t2 = getResult("Track",k);
			if(s1==s2 && t1==t2)
			{
				IJ.deleteRows(j,j); 
				area = a1+a2;
				setResult("Area",k,area);
				mean = (m1+m2)/2;
				setResult("Mean",k,mean);
				x = (x1+x2)/2;
				setResult("X",k,x);
				y = (y1+y2)/2;
				setResult("Y",k,y);
				d = (d1+d2)/2;
				setResult("Displacement",k,d);
				da = a/a1;
				setResult("Area Change",k,da);		
				updateResults();
				obsolete = newArray(k,j);
				roiManager("select",obsolete);
				index = getInfo("roi.name");
				roiManager("Combine");
				roiManager("Add");
				roiManager("select",roiManager("count")-1);
				roiManager("Rename",index);
				roiManager("select",obsolete);
				roiManager("delete");
				roiManager("sort");
				setBatchMode(false);
			}
		}
	}	
	setBatchMode("Hide");
}

function drawTracks()
{
	newImage("Tracks", "16-bit black", image_width, image_height, 1);
	trid = getImageID;	
	Array.getStatistics(tracks, min, max);
	for(i=1;i<=max;i++)
	{
		xcoord = newArray(0);
		ycoord = newArray(0);
		for(j=1;j<nr;j++)
		{
			if(tracks[j]==i)
			{
			 	xcoord = Array.concat(xcoord, getResult("X",j));
			 	ycoord = Array.concat(ycoord, getResult("Y",j));			 	
			}
		}
		selectImage(trid);
		makeSelection("polyline", xcoord, ycoord);
		setForegroundColor(i,i,i); 
		if(selectionType()!=-1)run("Draw");
	}
	run("Select None");
	resetMinAndMax();
	run("glasbey");
	return trid;
}

function colorCode(im_id, tracks)
{
	Array.getStatistics(tracks, min, max);	
	nr = nResults;
	rmc = roiManager("count");
	if(nr!=rmc)exit("Mismatch between Result nr and ROI nr");
	selectImage(im_id);
	run("Duplicate...", "title=Coded duplicate");
	coded_id = getImageID;
	selectImage(coded_id);
	Stack.getDimensions(image_width, image_height, channels, slices, frames);
	if(bitDepth() == 8)run("16-bit");
	Stack.setChannel(channels);
	run("Add Slice", "add=channel");
	Stack.setChannel(channels+1);
	for(i=0;i<rmc;i++)
	{
		code = tracks[i];
		roiManager("select",i);
		run("Set...", "value="+code+" slice");	 	
	}
	run("Select None");
	run("glasbey_on_dark");
	Stack.setDisplayMode("composite");
	for(c=1;c<=channels;c++)
	{
		Stack.setChannel(c);
		resetMinAndMax();
		run("Grays");
	}
	Stack.setChannel(channels+1);
	Stack.setFrame(frames);
	resetMinAndMax();
	toggleOverlay();
	return coded_id;
}

function nearestNeighbour(x,y,z,md)
{
	d = 100000; 
	bestfit = -1;
	nr = nResults;
	for(k=0;k<nr;k++)
	{
		if(getResult("Frame",k)==z)
		{
			x2 = getResult("X",k);
			y2 = getResult("Y",k);
			di = sqrt(pow(x-x2,2)+pow(y-y2,2));
			if(di<d && di<md)
			{
				d = di;
				bestfit = k; 
			}
		}
	}
	nn = newArray(bestfit,d);
	return nn;
}

function sample(id)
{
	// this function lists all relevant parameters  in  one results table, including the mean intensities of different channels
	selectImage(id);
	run("Properties...", " unit=pixel pixel_width=1 pixel_height=1");
	run("Set Measurements...", "area mean shape centroid stack redirect=None decimal=4");
	nr 			= nResults;
	tracks 		= newArray(nr);
	dists 		= newArray(nr);
	frame_nr 	= newArray(nr);
	area		= newArray(nr);
	means 		= newArray(0); //ccollect all means of all channels in one array (which will be of size nr * channels)
	area		= newArray(nr);
	// get track data from original results
	for(i=0;i<nr;i++)
	{
		tracks[i] 		= getResult("Track",i);
		dists[i] 		= getResult("Displacement",i);	
		frame_nr[i] 	= getResult("Frame",i);
	}
	// get means from channel 1
	for(c = 0; c < channels;  c++)
	{
		run("Clear Results");
		roiManager("Deselect");
		selectImage(id);
		Stack.setChannel(c+1);
		roiManager("Measure");
		for(i=0;i<nr;i++)
		{	
			means = Array.concat(means,getResult("Mean",i));
		}
	}
	// X, Y and Shape metrics should be retained from last measure command, other metrics are appended
	for(i=0;i<nResults;i++)
	{
		setResult("Track",i,tracks[i]);
		setResult("Displacement",i,dists[i]);	
		setResult("Frame",i,frame_nr[i]);
		for(c = 0; c < channels;  c++)
		{
			setResult("Mean_"+c+1,i,means[c*nr+i]);
		}
	}
	updateResults;
	return tracks;
}

function listTracks()
{
	if(isOpen("TrackData")){selectWindow("TrackData");run("Close");}
	nr 		= nResults;
	tracks 	= newArray(nr);
	times 	= newArray(nr);
	meana 	= newArray(nr);
	meanb 	= newArray(nr);
	// get track data from original results
	for(i=0;i<nr;i++)
	{
		tracks[i] 	= getResult("Track",i);
		times[i] 	= getResult("Frame",i);
		meana[i] 	= getResult("Mean_1",i);
		meanb[i] 	= getResult("Mean_2",i);
	}
	//	convert to track data
	nr = nResults;
	Array.getStatistics(tracks, min, max);
	Array.getStatistics(times, mint, maxt);
	// Create a new results table with metric organized per track
	IJ.renameResults("Results","Raw Data");
	// first fill matrix with zeros to avoid errors with jumping times
	for(i=0;i<nr;i++)
	{
		for(j=0;j<maxt;j++)
		{
			setResult("Track_"+tracks[i]+"_Mean_1",j,0);
			setResult("Track_"+tracks[i]+"_Mean_2",j,0);
		}
	}	
	updateResults;
	for(i=0;i<nr;i++)
	{
		setResult("Track_"+tracks[i]+"_Mean_1",times[i]-1,meana[i]);
		setResult("Track_"+tracks[i]+"_Mean_2",times[i]-1,meanb[i]);
	}
	updateResults;
	IJ.renameResults("Results","TrackData");
	IJ.renameResults("Raw Data","Results");
}

function add()
{
	if(selectionType==-1)waitForUser("Draw Region of Interest");
	setBatchMode(true);
	if(Stack.isHyperstack)Stack.getPosition(channel, slice, frame);
	else frame = getSliceNumber();
	getSelectionBounds(x, y, rw, rh);
	xc = x+rw/2;
	yc = y+rh/2;
	nr = nResults;
	track = getNumber("Track",0);
	roiManager("Add");
	roiManager("select",roiManager("count")-1);
	roi = Roi.getName;
	roiManager("Remove Channel Info");
	nf = d2s(frame,0);
	while(lengthOf(nf)<4)nf = "0"+nf; 
	nroi = ""+nf+substring(roi,4,14);
	roiManager("Rename",nroi);
	setResult("X",nr,xc); 
	setResult("Y",nr,yc); 
	setResult("Track",nr,track); 
	setResult("Frame",nr,frame); 
	updateResults();
	roiManager("deselect");
	run("Select None");
	 setBatchMode(false);
}

function update_rois(im_id)
{
	rois = roiManager("count");
	run("Select None");
	selectImage(im_id);
	for(r=0;r<rois;r++)
	{
		roiManager("select",r);
		roiManager("Update");
	}
	run("Select None");
	roiManager("Remove Channel Info");
}

function toggleOverlay()
{
	run("Select None"); roiManager("deselect");
	roiManager("Show All without labels");
	if(Overlay.size == 0)run("From ROI Manager");
	else run("Remove Overlay");
}
