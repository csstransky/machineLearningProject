This is a demo version (© Copyright CRCSI http://www.crcsi.com.au/ ) of the following automatic building detetcion algorithm published in

1. M. Awrangjeb, M. Ravanbakhsh and C. S. Fraser, “Automatic detection of residential buildings using LIDAR data and multispectral imagery,” ISPRS Journal of Photogrammetry and Remote Sensing, 65(5), 457-467, Sept. 2010.
2. M. Awrangjeb, C. Zhang and C. S. Fraser, "Building Detection in Complex Scenes Thorough Effective Separation of Buildings from Trees," Photogrammetric Engineering & Remote Sensing (PE&RS), vol. 78(7), 729-745, July 2012.
3. M. Awrangjeb, M. Ravanbakhsh and C. S. Fraser, “Automatic Building Detection Using LIDAR Data and Multispectral Imagery,” Digital Image Computing: Techniques and Applications (DICTA 2010), 1-3 Dec 2010, Sydney, Australia.
4. M. Awrangjeb, M. Ravanbakhsh and C. S. Fraser, “Building Detection from Multispectral Imagery and LIDAR Data Employing A Threshold-Free Evaluation System,” The ISPRS Commission III symposium on Photogrammetric Computer Vision and Image Analysis (PCV 2010), vol. XXXVIII, Part 3A, pp. 49-55, Paris, France, September 1-3, 2010.

To run the program you need to run buildingDetection.m in MATLAB environment. In order to change the input/sample data set or to put your own input data set update the file input.txt. The output will be available at finalCandidates.txt file where each line corresponds to four consecutive points of a detected rectangular building or building part. 
To run simpy use buildingDetection(input.txt,output.txt);

This MATLAB version is relatively unstable and due to limitation of MATLAB memory usage the software cannot handle large data sets. For a more stable version and test & evaluation purposes we recommend to get a trial version of the Barista software available at http://www.baristasoftware.com.au/ which has been implemented in Visual Studio 2005.

This software is protected by the copyright laws (© Copyright CRCSI http://www.crcsi.com.au/ ). Unauthorised use of the whole or part of the software is prohibited. 

For research only purposes the source code may be used. 

Note that the C/C++ source code and executable are available with the Barista Software at https://bitbucket.org/crcsi/baristasource/overview. This was developed by the University of Melbourne and CRC for Spatial Information and please contact info@crcsi.com.au for more information.
