#include <Adafruit_MPU6050.h>
#include <SD.h>
#include <Servo.h>



Adafruit_MPU6050 mpu; //initialize IMU

//declaring 2 servos for yaw and pitch angles, should change to a form of quaternion in the future
Servo servoX; 
Servo servoY;

//initialize yaw and pitch read by imu
double angleX = 0;
double angleY = 0;

double dt = 0.1; //time step in seconds

typedef struct PID{
    double kp;
    double ki;
    double kd;
    double sp;
    double error_last;
    double integral_error;
    double dt;

    PID(double k_p, double k_i, double k_d, double target, float time_step) {
        kp = k_p;
        ki = k_i;
        kd = k_d;
        sp = target;
        dt = time_step;
        error_last = 0;
        integral_error = 0;
    }
};


PID myPID(0.45,0.,0.1, 0, dt);

//initialize deflection to center the tvc
double servXmid = 40;  
double servYmid = 21; 

double servXlim = 20;
double servYlim = 20;

int b = 1;
int c = -1;
int ha = 1;

int commandCon = 2; // command int 1 is zero, 2 is center, 3 is dance, 4 is pid control

void setup() {
    Serial.begin(115200);

    Serial.setTimeout(0);

    while (!Serial)
        delay(10);

    if (!mpu.begin()) {
        Serial.println(F("Failed to find MPU6050 chip"));
    }

    mpu.setAccelerometerRange(MPU6050_RANGE_16_G);
    mpu.setGyroRange(MPU6050_RANGE_500_DEG);
    mpu.setFilterBandwidth(MPU6050_BAND_21_HZ);

    

    servoX.attach(5);  
    servoY.attach(8);

    center();

    delay(500);

    Serial.println(F("initialization done."));
}


void loop() {

    sensors_event_t a, g, temp;
    mpu.getEvent(&a, &g, &temp);


    angleX = angleX + rad2deg(g.gyro.x) * dt;
    angleY = angleY + rad2deg(g.gyro.y) * dt;

/*
    Serial.print("X: ");
    Serial.print(angleX); 
    Serial.print("  ");
    Serial.print("Y: ");
    Serial.println(angleY);
 */  
    switch(Serial.readString().toInt()){
        case 1:
            commandCon = 1;
        break; 
        case 2:
            commandCon = 2;
        break; 
        case 3:
            commandCon = 3;
        break; 
        case 4:
            angleX = 0;
            angleY = 0;
            center();
            commandCon = 4;
        break; 
    }

    switch(commandCon){
        case 1:
            zero();
        break; 
        case 2:
            center();
        break; 
        case 3:
            dance();
            delay(150);
        break; 
        case 4:
            control();
        break; 
    }
    
    delay(dt*500);
}

//centers servos
void center(){
    servoX.write(servXmid);
    servoY.write(servYmid);
}

//PID control based off IMU
void control(){
    deflectServX(getPID(&myPID, angleX, dt));
    deflectServY(getPID(&myPID, angleY, dt));
    //Serial.println(getPID(&myPID, angleX, dt) + servXmid);
    //Serial.println(getPID(&myPID, angleY, dt) + servYmid);
}

//servos to 0 degrees
void zero(){
    servoX.write(0);
    servoY.write(0);
}

//makes servos dance
void dance(){
    deflectServX(15 * b);
    deflectServY(15 * c);
    
    if(ha == 1){
        b = b * -1;
    }
    else{
        c = c * -1;
    }
    ha = ha * -1;
    
}

void deflectServX(double d){
    if(d > servXlim){
        servoX.write(servXmid + servXlim);
    }
    else if(d < -servXlim){
        servoX.write(servXmid - servXlim);
    }
    else{
        servoX.write(servXmid + d);
    }
}

void deflectServY(double d){
    if(d > servYlim){
        servoY.write(servYmid + servYlim);
    }
    else if(d < -servYlim){
        servoY.write(servYmid - servYlim);
    }
    else{
        servoY.write(servYmid + d);
    }
}

double getPID(PID *pid, double angle, double dt) {
    // IMPLEMENT PID.PY here with P = 0.5, I = 0.1, D = 0.5 and +/-15 degrees max and min
    double error = pid->sp - angle;
    double derivative_error = (error - pid->error_last) / dt;
    pid->integral_error += error * dt;
    double output = pid->kp * error + pid->ki * pid->integral_error + pid->kd * derivative_error;
    pid->error_last = error;

    return output;
}


double getServoOutput(double actuatorCom) {
    // returns the servo output given a certain angle deflection by getPID()
    return;
}

double rad2deg(double radian)
{
    double pi = 3.14159;
    return(radian * (180 / pi));
}