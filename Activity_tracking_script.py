import cv2
import datetime
import imutils
import numpy as np
import os
from collections import deque


##create two numpy arrays to place over images
#NB pixel size is written as 360x640 in np but would be 640x360 in cv2
imgL = np.zeros((360, 640, 3), np.uint8) #create a numpy array 360x640 = pixel size of video
imgR  = np.zeros((360, 640, 3), np.uint8) #create a numpy array 360x640 = pixel size of video

#open and read the txt file of video names
Batch = '3002' #tadpole "batch number" to process
myfile = open('Activity_videos.txt')#open this file
myfile = myfile.readlines() #read all the lines of this file
for line in myfile[1:]: #ignore header line
    line = line.strip() #removes unreadable characters from end of lines
    line = line.split() #ensures that lines are read line by line correctly
    columnBatchNo = 1 #specify the column number containing the batch numbers
    if line[columnBatchNo] != Batch: #if 'columnBatchNo' is not equal to the 'Batch' number specified then ignore this line and continue
        continue
    columnID = 4 #specify the column of tadpoleIDs
    columnTrial = 6 #specify  the column of trial number per tadpoleID
    columnvid = 12 #specify the column containing the video directory location
    columnLR = 9 #specify the column which states if tadpole on left 'L' or right 'R'
    #next few columns specify the XY coordinates of the tadpole containers
    columnLXtop = 16 #L tadpole top left X
    columnLYtop = 17 #L tadpole top left y
    columnLXbottom = 18 #L tadpole bottom right x
    columnLYbottom = 19 #L tadpole bottom right y
    columnRXtop = 20 #R tadpole top left X
    columnRYtop = 21 #R tadpole top left Y
    columnRXbottom = 22 #R tadpole bottom right X
    columnRYbottom = 23 #R tadpole bottom left Y
    vidx = (line[columnvid]) #object "vid" refers to the column stating a video file's location
    colLR = (line[columnLR]) #oject "colLR" refers to the column stating if tadpole is on the left or the right
    rectL = cv2.rectangle(imgL, (int(line[columnLXtop]), int(line[columnLYtop])), (int(line[columnLXbottom]), int(line[columnLYbottom])), (255, 255, 255), -1) #object rectL defines the regtangle coordinates, colour and thickness of the tadpole containers located on the left
    rectR = cv2.rectangle(imgR, (int(line[columnRXtop]), int(line[columnRYtop])), (int(line[columnRXbottom]), int(line[columnRYbottom])), (255, 255, 255), -1) #object rectR defines the regtangle coordinates, colour and thickness of the tadpole containers located on the right
    
    ##set up the conditions required for the .txt files where the results will be stored
    results_fn = ('_'.join(('Activity_videos/Cammy_results/ACT', line[columnID], line[columnTrial]))) #name file "file_location"/ACT_tadpoleID_trial according to cooresponding lines in inputted .txt file
    suffix = '.txt' #ensure that the file is of.txt type 
    results_txt = os.path.join(results_fn + suffix) #join the suffix onto the end of the results_fn file name
    myresults = open(results_txt, 'w') #open the results file.  "w" = Overwite any exsisting contents.
    printheader = print('TadpoleID', 'TrialID', 'FrameNo', 'Xcentroid', 'Ycentroid', 'DistX', 'DistY', 'PixelDist', 'CumulPixelDist', sep='\t', file=myresults) #print header lines in results file, separate each column by "\t"
    myresults.close() #close results file

    ##set up a series of lists to store values in for use later
    #list of points to follow detected object movement (contails) 
    pts = deque(maxlen=10) #maxlen = object location in previous x frames
    #list of points to calculate distance between consecutive cx and cy points
    cxpts = deque(maxlen=2) #maxlen = 2. Only hold current object X position and last object X position
    cypts = deque(maxlen=2) #maxlen = 2. Only hold current object Y position and last object Y position
    #set up list to calculate cumulative distance travelled:
    dist_travelled = deque() #hold all cumulative distance values calculated

    ##Define the video
    cap = cv2.VideoCapture(vidx) #object "cap".  Use cv2.videocapture to read the frames in a video
    #Get information on video frames - allows manipulation of custom number of frames to read
    framecount = int(cap.get(cv2.CAP_PROP_FRAME_COUNT)) # object "framecount". count the number of frames in a video and convert into an integer

    ##Define the background subtraction method to be used
    sub_history = 500 #keep at 500 - how many prior frames are used to determine stationary objects
    sub_varThreshold = 150 #values usually work between 90-150 - the pixel threshold value to be counted as a "moving tadpole" object
    subtractor = cv2.createBackgroundSubtractorMOG2(history=sub_history, varThreshold=sub_varThreshold, detectShadows=True) #object "subtractor". using the cv2.createBackgroundSubtractorMOG2 method

    ##set up a while loop to loop through each frame in the video file and perform all of the following actions
def getmask(self):
    while True:
        frame = cap.read() #for each frame, read from the video capture file
        frame = frame[1] #handle the frame from the video capture file
        #If video is finished - break from the loop, otherwise continue
        if frame is None:
            myresults.close() #close the file when the video is finished
            break
        centroids_per_frame = [] #create an empty list to store the number of centroids per frame in

        ##Define whether tracking will be for the left or right tadple
        #creates the Region of interest for a frame (ROI)
        if colLR == "R": 
           frame = cv2.bitwise_and(frame, rectR) #ROI = right tadpole
        if colLR == "L":
            frame = cv2.bitwise_and(frame, rectL) #ROI = left tadpole

        
        ##Define some attributs which will appear on final frame result
        text = 'NO TADPOLES' #default text to be displayed on video process window
        frame_pos = cap.get(cv2.CAP_PROP_POS_FRAMES) #for each frame, print the frame number

        ##Convert the frame to greyscale colour base
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        ##Blur the frame
        #kernal size indicates degree of blur
        #defalt kernal size should be 3,3
        blur_kernal = (5,5) 
        blur = cv2.GaussianBlur(gray, blur_kernal, 0) 
        
        #apply the results of the subtractor to the processed image
        mask = subtractor.apply(blur, learningRate=0.05) #learning rate = 0.05
        
        ##Construct a "mask" for the thresholded image
        minthresh = 60 #original = 60
        maxthresh = 255 #original = 255
        mask = cv2.threshold(mask, minthresh, maxthresh, cv2.THRESH_BINARY)[1] 
        #perform a series of dilations to remove any small blobs left in the mask
        mask = cv2.dilate(mask, None, iterations=3) #diltion - fills in any gaps in thresholded image #iterations = 3
        return mask, frame
        

    
        ##Place some information of the room status on frame
        #put tracking status and date + time on each frame
        #put info on frame images
        cv2.putText(getmask(frame), "Tracking status: {}".format(text), (10, frame.shape[0] - 25), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 255), 2)
        cv2.putText(frame, datetime.datetime.now().strftime("%A %d %B %Y %I:%M:%S%p"), (10, frame.shape[0] - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.35, (0, 0, 255), 1)
        #Put tadpole ID, frame number and cumulative distance travelled information on frame
        #put info on frame images
        cv2.rectangle(frame, (10, 2), (150,20), (255,255,255), -1) #rectangle location and colour
        cv2.putText(frame, "TadpoleID: {}".format(line[columnID]), (15, 15), cv2.FONT_HERSHEY_SIMPLEX, 0.5 , (0,0,0)) #retrieve frame number from cap and place in rectangle
        cv2.rectangle(frame, (10, 40), (150,60), (255,255,255), -1) #rectangle location and colour
        cv2.putText(frame, "Frame No: {}".format(int(frame_pos)), (15, 55), cv2.FONT_HERSHEY_SIMPLEX, 0.5 , (0,0,0)) #retrieve frame number from cap and place in rectangle
        cv2.rectangle(frame, (10, 80), (150,100), (255,255,255), -1) #rectangle location and colour for frame number
        #cv2.putText(frame, "Dist (pixels): {}".format(cumul_dist_travelled), (15, 95), cv2.FONT_HERSHEY_SIMPLEX, 0.3 , (0,0,0)) #retrieve frame number from cap and place in rectangle
    

        ##decide what videos to show

        #single videos
        #cv2.imshow('frame', frame) #video containing only specified ROI on the original frame
        #cv2.imshow('blur', blur) #video containing only blurred image
        #cv2.imshow('mask', mask) #video containing object mask

        #multiple videos 
        #two show multiple videos, both videos must be displayed in the same channel
        #videos in greyscale must be converted to colour to be displayed next to frame video
        #convert greyscale to colour
        blur_2_colour = cv2.cvtColor(blur, cv2.COLOR_GRAY2BGR)
        mask_2_colour = cv2.cvtColor(mask, cv2.COLOR_GRAY2BGR)
        #stack the images together
        hstack_frame_blur = np.hstack((frame, blur_2_colour))
        hstack_frame_mask = np.hstack((frame, mask_2_colour))

        #cv2.imshow('frame_blur', hstack_frame_blur)
        cv2.imshow('frame_mask', hstack_frame_mask)
        
    
                ##how long (in millisecs) does each frame appear for?
        key = cv2.waitKey(60) & 0xFF

        ##To skip a video press:
        if key == ord('q'):
            break
 
    #close files and windows after analysis is complete        
    cap.release() #release the video
    cv2.destroyAllWindows() #destroy any video windows open        
    myresults.close() #close the file when a video is finished or no more videos are playing


