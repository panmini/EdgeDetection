#include "BMP.h"

BMP::BMP()
{
    imageinfo = new ImageInfo;
    fileinfo = new FileInfo;
}
void BMP::readImage(string filename)
{
    ifstream file(filename , ios::binary); //opening file

    char filetemp[14];                     //store file informations
    char imagetemp[40];

    file.read(filetemp , sizeof(filetemp)); //read first 14 byte
    fileinfo->setFileHeader(filetemp);      //send allData
    file.read(imagetemp , sizeof(imagetemp));
    imageinfo->setImageHeader(imagetemp);
    fileinfo  ->  setOffSet(54);

    file.seekg(54);    //move cursor offset

    padding = 0;

    while((imageinfo->getWidth()*3+padding)%4 != 0) padding++; //calculate padding

    data = new BYTE[imageinfo -> getBiSize()];  //determine data size
    BYTE * pointerOfData = data;                //point data
    BYTE buffer[imageinfo->getWidth()*3];       //temprature memory

    for(int i=0;i<(int)(imageinfo->getBiSize()/sizeof(buffer));i++)
	{
        file.read((char*)buffer,sizeof(buffer));        //read first row image
        memcpy(pointerOfData,buffer,sizeof(buffer));    //copy buffer adresses
        pointerOfData += sizeof(buffer);                //move pointer
        file.read((char*)buffer,padding);               //move cursor
    }

    file.close();
}

void BMP::saveGrayScale(string filename)
{
    grayData = new BYTE[imageinfo->getBiSize()/3];
    BYTE * iterator = grayData;
    BYTE* pointerOfData = grayData;


    for(int i=0;i<(int)imageinfo->getBiSize();i+=3)
	{
        *iterator = BYTE((data[i]*0.21+data[i+1]*0.72+data[i+2]*0.07));
        iterator++;
    }

    ofstream file(filename + ".bmp" , ios::binary);
    file.write((char *)fileinfo->getAllHeader() , 14);
    file.write((char *)imageinfo->getAllInfo() , 40);


    unsigned int pixelNumber = 0;

    for(int i=0;i<(int)imageinfo->getBiSize()/3;i++)
	{
        file.write((char*)pointerOfData,1);
        file.write((char*)pointerOfData,1);
        file.write((char*)pointerOfData++,1);

        pixelNumber++;

        if(pixelNumber == imageinfo->getWidth())
		{
            BYTE pad = 0;
            for(int i=0;i<padding;i++) file.write((char*)&pad,1);
            pixelNumber = 0;
        }
    }

    file.close();
}
void BMP::saveGrayScale(BYTE * temp , string filename)
{
    BYTE * pointerOfData = temp;

    ofstream file(filename + ".bmp" , ios::binary);
    file.write((char *)fileinfo->getAllHeader() , 14);
    file.write((char *)imageinfo->getAllInfo() , 40);

    unsigned int pixelNumber = 0;

    for(int i=0;i<(int)imageinfo->getBiSize()/3;i++)
	{
        file.write((char*)pointerOfData,1);
        file.write((char*)pointerOfData,1);
        file.write((char*)pointerOfData++,1);

        pixelNumber++;

        if(pixelNumber == imageinfo->getWidth())
		{
            BYTE pad = 0;
            for(int i=0;i<padding;i++) file.write((char*)&pad,1);
            pixelNumber = 0;
        }
    }
    file.close();
}

void BMP::saveImage(string filename)
{
    ofstream file(filename + ".bmp" , ios::binary);
    file.write((char *)fileinfo->getAllHeader() , 14);
    file.write((char *)imageinfo->getAllInfo() , 40);

    BYTE * pointerOfData = data;
    unsigned int pixelNumber = 0;

    for(int i=0;i<(int)imageinfo->getBiSize();i++)
	{
        file.write((char*)pointerOfData++,1);

        pixelNumber++;
        if(pixelNumber >= imageinfo->getWidth() * 3)
		{
            BYTE pad = 0;
            for(int i=0;i<padding;i++) file.write((char*)&pad,1);
            pixelNumber = 0;
            i += padding;
        }
    }

    file.close();
}

BYTE* BMP::Conv(BYTE* raw, int width, int height, float* convMatrix, int cW) //-> çerçeveyi sildim
{
	int w =  width + cW - 1;
	int h =  height + cW - 1;

	BYTE* zero = new BYTE[(w) * (h)];

	for (int i = 0; i < w * h; i++)
		zero[i] = 0;

	for(int i = 0; i< cW/2 ;i++)
		for (int j = 0; j < w ; j++)
		{
			zero[i*(w) + j] = 0;
			zero[((h) - i - 1)*(w) + j] = 0;
		}

	for (int i = 0; i < (h); i++)
		for (int j = 0; j < cW/2; j++)
		{
			zero[i*(w) + j] = 0;
			zero[i*(w) + j+ (w) - cW/2] = 0;
		}

	for (int i = cW/2  ; i < (h)-cW/2 ; i++)

		for (int j = cW/2 ; j < (w)-cW/2 ; j++)
		{
			zero[i * (w) + j] = *raw;
			raw++;
		}

	float norm = 0;
	for (int i = 0; i < cW*cW; i++)
		norm += convMatrix[i];

	if ((int)norm == 0) norm = 1;
	int sum;

	BYTE* borderless = new BYTE[width*height];
	for (int i = 0; i < width * height; i++)
		borderless[i] = 0;

	int a = 0;
	for (int i = 0; i <= h - cW; i++)
		for (int j = 0; j <= w - cW; j++)
		{
			sum = 0;
			for (int k = 0; k < cW; k++)
				for (int m = 0; m < cW; m++)
					sum += int(zero[(w * k + m) + j + (i*(w))] * convMatrix[cW * k + m]);
			borderless[a] = BYTE((int)abs(sum / norm));	a++;
		}
	delete[] zero;

	return borderless;
}

BYTE* BMP::SobelFiltering(BYTE* img, int w,int h, float* edge_direction)
{

    float Gx[9]={-1.0,0.0,1.0,
                -2.0,0.0,2.0,
                -1.0,0.0,1.0};

    float Gy[9]={1.0 ,2.0 ,1.0,
                0.0  ,0.0  ,0.0,
                -1.0  ,-2.0  ,-1.0};

    BYTE* sobelX=Conv(img, w,h,Gx,3);       //****************
    BYTE* sobelY=Conv(img, w,h,Gy,3);
    saveGrayScale(sobelX,"sobelx");
	saveGrayScale(sobelY,"sobely");


    float angle;
    float *edgeMagnitude=new float[w*h];
    
    float max=0.0;
    BYTE* temp=new BYTE[w*h];
    for(int i=0; i< w;i++)
        for(int j=0;j<h;j++)
        {
            angle=0.0;
            edgeMagnitude[i+j*w]=sqrt(sobelX[i+j*w]*sobelX[i+j*w]+sobelY[i+j*w]*sobelY[i+j*w]);
            max = edgeMagnitude[j * w + i] > max ? edgeMagnitude[j * w + i] : max;
            ///////angle calc
            if ((sobelX[i+j*w] != 0.0) || (sobelY[i+j*w] != 0.0))
            {
                angle = atan2(sobelY[i+j*w], sobelX[i+j*w]) * 180.0 / 3.1415926535897; // 180 / pi rad to degree
                
            }
            else
            {
                angle = 0.0;
            }

            if (((angle > -22.5) && (angle <= 22.5)) || ((angle > 157.5) && (angle <= -157.5)))
            {
                edge_direction[j * w + i] = 0;
            }
            else if (((angle > 22.5) && (angle <= 67.5)) || ((angle > -157.5) && (angle <= -112.5)))
            {
                edge_direction[j * w + i] = 45;
            }
            else if (((angle > 67.5) && (angle <= 112.5)) || ((angle > -112.5) && (angle <= -67.5)))
            {
                edge_direction[j * w + i] = 90;
            }
            else if (((angle > 112.5) && (angle <= 157.5)) || ((angle > -67.5) && (angle <= -22.5)))
            {
                edge_direction[j * w + i] = 135;
            }

            ///////
        }

        for (int x = 0; x < h; x++)
        {
            for (int y = 0; y < w; y++)
            {
                edgeMagnitude[x * w + y] = 255.0f * edgeMagnitude[x * w + y] / max;
                temp[x*w+y]=edgeMagnitude[x*w+y];
            }
        }

    NonMaxSuppression(edge_direction,edgeMagnitude,temp,w,h);

	double lowThresholdRatio=0.5, highThresholdRatio=0.183;


	BYTE lowThreshold, highThreshold;


	highThreshold	= max*highThresholdRatio;
	lowThreshold	= highThreshold*lowThresholdRatio;

	cout <<"low"<<(int)lowThreshold<<"\t high"<<(int)highThreshold<<"\n\n";

	Hysteresis(lowThreshold,highThreshold,w,h,temp);

    return temp;
}

void BMP::NonMaxSuppression(float* edge_direction, float* edge_magnitude,BYTE* value, int width, int height)
{
    saveGrayScale(value,"sobel");
	float pixel_1 = 0;
	float pixel_2 = 0;
	float pixel;

	for (int x = 1; x < height - 1; x++) {
		for (int y = 1; y < width - 1; y++) {
			if ((edge_direction[x * width + y] >= 0 && edge_direction[x * width + y] < 22.5)        //***************** == 0 / 45 / 90 / 135 yazacasın
			|| (edge_direction[x * width + y] >= 157.5 && edge_direction[x * width + y] <=180))
			{
				pixel_1 = edge_magnitude[(x ) * width + y+1];
				pixel_2 = edge_magnitude[(x ) * width + y-1];
			}
			else if (edge_direction[x * width + y] >= 22.5 && edge_direction[x * width + y] < 67.5) {
				pixel_1 = edge_magnitude[(x + 1) * width + y - 1];
				pixel_2 = edge_magnitude[(x - 1) * width + y + 1];
			}
			else if (edge_direction[x * width + y] >= 67.5 && edge_direction[x * width + y] < 112.5) {
				pixel_1 = edge_magnitude[(x+1) * width + y ];
				pixel_2 = edge_magnitude[(x-1) * width + y];
			}
			else if (edge_direction[x * width + y] >= 112.5 && edge_direction[x * width + y] < 157.5) {
				pixel_1 = edge_magnitude[(x - 1) * width + y - 1];
				pixel_2 = edge_magnitude[(x + 1) * width + y + 1];
			}
			pixel = edge_magnitude[x * width + y];
			if ((pixel >= pixel_1) && (pixel >= pixel_2)) {
                value[x*width+y]=pixel;

			} else {

                value[x*width+y]=0;
			}
		}
	}
    saveGrayScale(value,"nonmax");    

}

void BMP::Hysteresis(BYTE lowThreshold, BYTE highThreshold, int w, int h, BYTE* value)
{
	for(int i=0; i<w;i++)
		for(int j=0;j<h;j++)
		{
			if(value[j*w+i] >= highThreshold) 		value[j*w+i]=255;
			else if(value[j*w+i] < lowThreshold)	value[j*w+i]=0;
			//else value[j*w+i]=128;
		}
	saveGrayScale(value,"thres");

	for(int i=1; i<w-1;i++)
		for(int j=1;j<h-1;j++)
		{
            if( value[i+j*w] != 0 )
			    if( value [i-1+(j-1)*w] == 255 || value [i+(j-1)*w] == 255 || value [i+1+(j-1)*w] == 255
			    || value [i-1+j*w] == 255 || value [i+1+j*w] == 255
			    || value [i-1+(j+1)*w] == 255 || value [i+(j+1)*w] == 255 || value [i+1+(j+1)*w] == 255
			    )
			       	value[i+j*w]=255;
			    else value[i+j*w]=0;
    	}

	saveGrayScale(value,"Hysteresis");
}

void BMP::houghTransform(BYTE* value, int w, int h, float* edge_direction)
{   /////Set Accu values
    double hough_h = (sqrt( (double)(h*h + w*w)));      //*****hough_h sil
    int accu_h = hough_h ;
    int accu_w = 180;

    double DEG2RAD = 3.14159265/180.0;
    BYTE* accu = new BYTE[accu_h*accu_w];
    double d;

    for(int i=0;i<accu_h*180;i++)
        accu[i]=0;

    for(int y= 0;y<h;y++)
    {
        for(int x= 0;x<w;x++)
        {
            if( ( value [ (y*w) + x ] ) == 255 )
            for (double t = 0; t < 180; t++)
            {
		        double theta = ( t ) * 3.14159265/180.0;
                d=((double)x)* cos( theta ) + (double)y* sin( theta );
                    
                accu[ (int)(round(d)*accu_w + t) ]++;
            }            
        }
    }

    //for saving accu set image info to accus values
        imageinfo->setHeight(accu_h);
        imageinfo->setWidth(180);
        
        padding = 0;
        while(((imageinfo->getWidth())*3+padding)%4 != 0) padding++;

        DWORD imageSize = accu_h * accu_w * 3;
        fileinfo  ->setFileHeader((char*)fileinfo->getAllHeader());
        fileinfo  ->setSize(imageSize);
        imageinfo -> setImageHeader((char *)imageinfo->getAllInfo());
        imageinfo -> setSize(imageSize - 54);

        saveGrayScale(accu, "accu"); 
    /////////
    imageinfo->setHeight(h);
    imageinfo->setWidth(w);
        
    padding = 0;
    while(((imageinfo->getWidth())*3+padding)%4 != 0) padding++;

    imageSize = h * w * 3;
    fileinfo  ->setFileHeader((char*)fileinfo->getAllHeader());
    fileinfo  ->setSize(imageSize);
    imageinfo -> setImageHeader((char *)imageinfo->getAllInfo());
    imageinfo -> setSize(imageSize - 54);

    BYTE* hough=new BYTE[w*h];
    for(int i=0;i<w*h;i++)
    {
        hough[i]=grayData[i];
        grayData[i]=255;
    }

    int max;
    int iter;
    int x1, y1, x2, y2;
    x1 = y1 = x2 = y2 = 0;

    for(int i=0;i< 180 ;i++)    //*********************
    {
        max=0;
        for(int k=0;k< 180 ;k++)
            for(int j=0;j<accu_h;j++)
                if(accu[j*accu_w+k] > max)
                {
                    max = accu[j*accu_w+k];
                }

        for(int j=0;j<accu_h;j++)
            if(accu[j*accu_w+i] > max * 0.4)
            {
                iter = j;           
            
                if ( i == 45 || i == 135 || i == 0 ||  i == 90)
                {
        
                    if(i >=45  && i <= 135)
		            {
        	            //y = (r - x cos(t)) / sin(t)
			            x1 = 0;
			            y1 = ((double)(iter) - ((x1  ) * cos((i ) * DEG2RAD ))) / sin((i ) * DEG2RAD);
			            x2 = w - 0;
			            y2 = ((double)(iter) - ((x2 ) * cos((i ) * DEG2RAD))) / sin((i ) * DEG2RAD);
            
		            }
		            else 
		            {
		    	        //x = (r - y sin(t)) / cos(t);
		    	        y1 = 0;
                        x1 = ((double)(iter) - ((y1  ) * sin((i ) * DEG2RAD))) / cos((i ) * DEG2RAD)  ;
		    	        y2 = h - 0;
			            x2 = ((double)(iter) - ((y2 ) * sin((i ) * DEG2RAD))) / cos((i ) * DEG2RAD) ;
                    }
                }       
      
            Line(x1,y1,x2,y2);
            }
    }

    saveGrayScale(grayData,"houghTransform");
    grayData=hough;   
}

void BMP::Line( float x1, float y1, float x2, float y2)
{   
    cout <<"x1 -- y1 -- x2 -- y2\n"<<x1<<"\t"<<y1<<"\t"<<x2<<"\t"<<y2<<"\n";
    int w=imageinfo->getWidth();
    int h=imageinfo->getHeight();
    // Bresenham's line algorithm
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if(steep)
    {
        swap(x1, y1);
        swap(x2, y2);
    } 

    if(x1 > x2)
    {
        swap(x1, x2);
        swap(y1, y2);
    }

    const float dx = x2 - x1;
    const float dy = fabs(y2 - y1);

    float error = dx / 2.0f;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for(int x=(int)x1; x<maxX; x++)
    {
        if(steep)    
            grayData[x*w+y]=0;    
        else    
            grayData[x+w*y]=0;    

        error -= dy;
        if(error < 0)
        {
            y += ystep;
            error += dx;
        }
    }
}

void BMP::finalisation()
{
	float gaussian3[9]    = {   0.25,0.5,0.25,
                                0.5,1,0.5,
                                0.25,0.5,0.25};

    float convMatrix5[25] = {   0.0625,	0.125,	0.25,	0.125,	0.0625,
								0.125,	0.25,	0.5,	0.25,	0.125,
								0.25,	0.5,	1,		0.5,	0.25,
								0.125,	0.25,	0.5,	0.25,	0.125,
								0.0625,	0.125,	0.25,	0.125,	0.0625 }; //7

	float convMatrix7[49] = {	0.015625 ,	0.03125,	0.0625,	0.125,	0.0625,	0.03125,0.015625,
								0.03125,	0.0625,		0.125,	0.25,	0.125,	0.0625,	0.03125,
								0.625,		0.125,		0.25,	0.5,	0.25,	0.125,	0.625,
								0.125,		0.25,		0.5,	1,		0.5,	0.25,	0.125,
								0.625,		0.125,		0.25,	0.5,	0.25,	0.125,	0.625,
								0.03125,	0.0625,		0.125,	0.25,	0.125,	0.0625,	0.03125,
								0.015625,	0.03125,	0.0625, 0.125,	0.0625,	0.03125,0.015625 };
                                
    float *edge_direction=new float[imageinfo->getWidth()*imageinfo->getHeight()];

    BYTE* gaussian=Conv(grayData, imageinfo->getWidth(), imageinfo->getHeight(), gaussian3,3);
    
    saveGrayScale(gaussian,"blur");

    BYTE* toHough=SobelFiltering(gaussian,imageinfo->getWidth(),imageinfo->getHeight(),edge_direction);

    houghTransform(toHough,imageinfo->getWidth(),imageinfo->getHeight(), edge_direction);

}
