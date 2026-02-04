import numpy as np
import time
import PID as pid
import Graphics as gh
import Physics as phys
import rocketConfig
import pandas as pd
#Units are SI units... m/s/radians/kg

time_ret = []
angles_ret = []
vert_pos_ret = []
hori_pos_ret = []

ori_labels = {"xlab" : "time(s)", "ylab" : "angle(theta)", "title" : "Theta"}
translate_labels = {"xlab" : "time(s)", "ylab" : "Pos Z", "title" : "Position Z"}
vertical_labels = {"xlab" : "time(s)", "ylab" : "Pos X", "title" : "Position X"}

graphs_dict = [ori_labels,translate_labels,vertical_labels]

graphics = gh.GraphicHandler()
graphics.graphsHandler(3,graphs_dict) #number of graphs and their labels

rocket = graphics.createAgent('black')
targe = graphics.createAgent('red')

vehicle = rocketConfig.vehicleProperties(0.25,0.1,0.5,np.deg2rad(15),np.deg2rad(150))
#mass(kg),mmoi(kg m^2),com2TVC(meters),servo_lim(rad),servo_rate_lim(rad/s)

#THRUST = 12 # Newtons


#initial state vector
state_vector = {"ax" : 0 ,"vx" : 0,"px" : 0,"az" : 0 ,"vz" : 1,"pz" : 0 ,"alpha" : 0.0,"omega" : 0.05,"theta" : 0.5}


rocket_phys = phys.threeDofPhysics(state_vector,vehicle.mass,vehicle.mmoi)

#Controller setup
controller = pid.PID(0.5,0.05,0.1,0)
#controller = pid.PID(.85,0.2,0.1,0)

controller.setLims(-10,10)#output limits
#our TVC is also limited by SERVO_LIMIT but we might want to change the the two independently

sim_time = 0.0
time_lim = 15.0
delta_t = 0.1
tvc_input = 0


def getThrust (time):
    datafile = open("PYTHON SIM(old)/motor.txt","r")
    lines = datafile.readlines()
    coun = 0
    time1 = 0
    thrust1 = 0
    time2 = 0
    thrust2 = 0
    while coun<len(lines) and time >= float(lines[coun].split(" ")[0]): 
        coun+=1
    if (coun == 0 or coun == len(lines)):
        return 0
    else :
        time1 = float(lines[coun-1].split(" ")[0])
        thrust1 = float(lines[coun-1].split(" ")[1].split("\n")[0])
        time2 = float(lines[coun].split(" ")[0])
        thrust2 = float(lines[coun].split(" ")[1].split("\n")[0])
        if (time == time1):
            return thrust1
        else:
            return (((thrust2-thrust1)/(time2-time1))*(time-time1)) + thrust1


if __name__ == "__main__":
	while(sim_time < time_lim):

		#Physics
		forces = rocket_phys.tvcPhysics(tvc_input,getThrust(sim_time),vehicle,delta_t) # Compute Forces	from TVC	
		rocket_phys.inputForces(forces,delta_t) # Apply Forces to body	
		graphics.moveAgent(rocket,rocket_phys.state_vector["pz"],rocket_phys.state_vector["px"])# Update Graphics and Plots
		graphics.rotateAgent(rocket,rocket_phys.state_vector["theta"])

		#Saving data we want to plot
		angles_ret.append(np.rad2deg(rocket_phys.state_vector['theta']))
		hori_pos_ret.append(rocket_phys.state_vector["pz"])
		vert_pos_ret.append(rocket_phys.state_vector["px"]) 
		graphs = [(time_ret,angles_ret),(time_ret,hori_pos_ret),(time_ret,vert_pos_ret)]

		#Real-time plotting is slow right now -> need to look into blit for matplotlib
		#graphics.updateGraphs(graphs)

		# Compute Control Logic
		tvc_input = controller.compute(rocket_phys.state_vector['theta'],delta_t)


		time_ret.append(sim_time)
		sim_time += delta_t

	graphics.showGraphs(graphs)

	#edd
