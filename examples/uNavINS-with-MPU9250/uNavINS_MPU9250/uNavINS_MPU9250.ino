/*
uNavINS_MPU9250.ino
Brian R Taylor
brian.taylor@bolderflight.com

Copyright (c) 2017 Bolder Flight Systems

Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
and associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "uNavINS.h"
#include "MPU9250.h"
#include "UBLOX.h"

// an MPU-9250 object on SPI bus 0 with chip select 24
MPU9250 Imu(SPI,24);
int status;
// a flag for when the MPU-9250 has new data
volatile int newData;
// a uNavINS object
uNavINS Filter;
UBLOX gps(4);
bool newGpsData;
unsigned long prevTOW;
gpsData uBloxData;
// timers to measure performance
unsigned long tstart, tstop;

void setup() {
  // serial to display data
  Serial.begin(115200);
  while(!Serial) {}

  gps.begin(115200);

  // start communication with IMU 
  status = Imu.begin();
  if (status < 0) {
    Serial.println("IMU initialization unsuccessful");
    Serial.println("Check IMU wiring or try cycling power");
    Serial.print("Status: ");
    Serial.println(status);
    while(1) {}
  }
  // setting SRD to 9 for a 100 Hz update rate
  Imu.setSrd(9);
  // enabling the data ready interrupt
  Imu.enableDataReadyInterrupt();
  // attaching the interrupt to microcontroller pin 1
  pinMode(27,INPUT);
  attachInterrupt(27,runFilter,RISING);
  Serial.println("Starting...");
}

void loop() {
  gps.read(&uBloxData);
  if (uBloxData.numSV > 5) {
    if (newData == 1) {
      newData = 0;
      tstart = micros();
      // read the sensor
      Imu.readSensor();
      // update the filter
      Filter.update(uBloxData.iTOW,uBloxData.velN,uBloxData.velE,uBloxData.velD,uBloxData.lat*PI/180.0f,uBloxData.lon*PI/180.0f,uBloxData.hMSL,Imu.getGyroY_rads(),-1*Imu.getGyroX_rads(),Imu.getGyroZ_rads(),Imu.getAccelY_mss(),-1*Imu.getAccelX_mss(),Imu.getAccelZ_mss(),Imu.getMagX_uT(),Imu.getMagY_uT(),Imu.getMagZ_uT());
      Serial.print(Filter.getPitch_rad()*180.0f/PI);
      Serial.print("\t");
      Serial.print(Filter.getRoll_rad()*180.0f/PI);
      Serial.print("\t");
      Serial.print(Filter.getYaw_rad()*180.0f/PI);
      Serial.print("\n");
      tstop = micros();
    }
  }
}

void runFilter() {
  newData = 1;
}
