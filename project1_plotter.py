import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    file_object = open(filename, 'r')
    lines = file_object.readlines()
    file_object.close()

    return lines

def plotter_solution(filenames): 

    for filename in filenames:   
        lines = read_file(filename)
        x = np.zeros(len(lines)) #x-values
        v = np.zeros(len(lines)) #will contain v-vector from file
        v_closed_form = np.zeros(len(lines)) #closed form solution for v. 

        i = 0
        for line in lines:
            x[i] = line.split()[0]
            v[i] = line.split()[1]
            v_closed_form[i] = line.split()[2]
            i +=1

        picturename = "blabla_%d.png" %(len(lines))
        plt.plot(x, v, label = "numerically calculated solution")
        plt.plot(x, v_closed_form, label = "closed-form solution")
        plt.legend()
        plt.title("n = %g" %len(x))
        plt.xlabel("grid points")  
        plt.ylabel("v")   
        plt.savefig(picturename)
        plt.show()

plotter_solution(["ownsolver_n_10.dat", "ownsolver_n_100.dat", "ownsolver_n_1000.dat"])

def relative_error(filename): 
    lines = read_file(filename)

    n = np.zeros(len(lines))
    error = np.zeros(len(lines))

    i = 0
    for line in lines: 
        n[i] = line.split()[0]
        error[i] = line.split()[1]
        i+=1

    plt.loglog(n, error, "-*", label = "relative error")
    plt.legend()
    plt.xlabel("log10(n)")
    plt.ylabel("log10(relative error)")
    plt.savefig("relative_error.png")
    plt.show()

relative_error("relative_error.dat")

    
    
