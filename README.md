# MODIS-SST-Granule-to-Lat-Lon

MATLAB script to reproject MODIS Sea Surface Temperature (SST) Level 2 data freely available through Oceancolor (https://oceancolor.gsfc.nasa.gov/). 

The Level 2 data is provided in 5-minute segements, called "granules". A MODIS granule covers an area of 2030 km along path by 2330 km across path.  The MATLAB script reprojects the granular data into a lat/lon projection using a nearest neighbor algorithm. 

If you find this script useful and use it in your work, please cite the following citation:

R. Cassotto, M. Fahnestock, J. M. Amundson, M. Truffer, and I. Joughin, “Seasonal and interannual variations in ice mélange and its impact on terminus stability, Jakobshavn Isbræ, Greenland,” Journal of Glaciology, vol. 61, no. 225, pp. 76–88, 2015, doi: 10.3189/2015jog13j235.

![image](https://github.com/user-attachments/assets/e18d8987-d3a3-4492-9c5c-cc1900aecc9d)
