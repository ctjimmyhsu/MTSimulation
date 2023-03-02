from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import random

##################################################

#savedir = 'C:/Users/brouhardlab/Desktop/'

def simulate(tottime, interval, vg, vs, numtubes, steps, rate, loc, nuctime, bypassflag):
    
    interval = float(interval)
    vg = float(vg)
    vs = float(vs)
    
    scale = 1/rate

    lengths = [0]*numtubes #Keep track of lengths
    
    states = [0]*numtubes # -1 = shrinking, 0 = nucleating, 1 = growing
    
    catdecisions = [0]*numtubes #Keep track of descision to catastrophe
    
    nucdecisions = [0]*numtubes #Keep track of descision to nucleate

    nucleationtimesteps = [0]*numtubes #To keep track of how long a MT has been trying to nucleate
    
    growthtimesteps = [0]*numtubes #To keep track of how long a MT has been growing

    lengthmt1 = [] #An example microtubule for output at the end
    
    probrange = np.arange(0,3600,1) #range to build probability distributions
    
    ##############################
    # CATASTROPHE PROBABILITIES
    
    ##
    gammadis = stats.gamma.cdf(probrange, steps, loc, scale)
    ##
    
    fig, ax = plt.subplots(1,1, figsize=(4, 3.5))
    ax = plt.subplot(1,1,1)
    
    plt.plot(probrange, gammadis, linewidth = 2, color = 'k')
    
    plt.grid(False)
    ax.set_facecolor('white')
    ax.tick_params(axis='both', color = 'k', length = 4, width = 1.5)
    plt.xlim(0,900)
    plt.ylim(0,1)
    plt.title('Lifetime IN', fontsize = 15, color = 'k')
    plt.xlabel('Lifetime (sec)', fontsize = 15, color = 'k')
    plt.ylabel('Cumulative Probability', fontsize = 15, color = 'k')
    plt.xticks(fontsize = 15, color = 'k')
    plt.yticks(fontsize = 15, color = 'k')
    plt.axhline(y = 0, color='k', linewidth=2.5)
    plt.axvline(x = 0, color='k', linewidth=2.5)
    #ax.text(0.1, 140, str(i*interval) + " sec", verticalalignment='center', horizontalalignment='left', fontsize = 15, color = "black")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #plt.savefig(savedir+'gammadis.jpg')
    plt.show()
    plt.tight_layout()
    
    ##############################
    # NUCLEATION PROBABILITIES
    
    if bypassflag != 1:
        
        ##
        nucldis = stats.expon.cdf(probrange, 0, nuctime)
        ##
    
        fig, ax = plt.subplots(1,1, figsize=(4, 3.5))
        ax = plt.subplot(1,1,1)

        plt.plot(probrange, nucldis, linewidth = 2, color = 'k')

        plt.grid(False)
        ax.set_facecolor('white')
        ax.tick_params(axis='both', color = 'k', length = 4, width = 1.5)
        plt.xlim(0,900)
        plt.ylim(0,1)
        plt.title('Nucleation IN', fontsize = 15, color = 'k')
        plt.xlabel('Time (sec)', fontsize = 15, color = 'k')
        plt.ylabel('Cumulative Probability', fontsize = 15, color = 'k')
        plt.xticks(fontsize = 15, color = 'k')
        plt.yticks(fontsize = 15, color = 'k')
        plt.axhline(y = 0, color='k', linewidth=2.5)
        plt.axvline(x = 0, color='k', linewidth=2.5)
        #ax.text(0.1, 140, str(i*interval) + " sec", verticalalignment='center', horizontalalignment='left', fontsize = 15, color = "black")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #plt.savefig(savedir+'nucldis.jpg')
        plt.show()
        plt.tight_layout()
    
    #########################################
    #########################################
    #  SIMULATION
    
    avglengths_lst = []
    lifetimes = []
    
    ### Increment by one time interval
    for i in np.arange(0,tottime,interval):
    ###
    
        avglength = 0 #Reinitialized for each timestep

        ### Simulate one microtubule at a time
        for m in np.arange(0,numtubes):
        ###

            ##########################
            #If reached starting point, decide when to nucleate

            if states[m] == 0:
                
                if bypassflag == 1:
                    
                    nucdecisions[m] = 0
                    
                else:
                
                    nucdecisions[m] = np.random.exponential(nuctime)

                states[m] = 1

            ##########################
            #If still waiting to nucleate, just increment
                    
            if states[m] == 1:
                
                nucleationtimesteps[m] += 1
                
                currnucleationtime = nucleationtimesteps[m]*interval
                
                if currnucleationtime > nucdecisions[m]:
                    
                    states[m] = 2
                    
                    nucleationtimesteps[m] = 0

            ##########################
            #If decided to nucleate, decide when to catastrophe

            if states[m] == 2:
                
                catdecisions[m] = np.random.gamma(steps,scale)

                states[m] = 3
                
            ##########################
            #Grow and check if you need to catastrophe

            if states[m] == 3:
                
                lengths[m] = lengths[m] + (vg/60)*interval
                
                growthtimesteps[m] += 1
                
                currlifetime = growthtimesteps[m]*int(interval)

                if currlifetime > catdecisions[m]:
                
                    lifetimes.append(currlifetime)
                    
                    states[m] = -1

                    growthtimesteps[m] = 0

            ##########################
            #If catastrophed, shrink until the end

            if states[m] == -1:
                
                lengths[m] = lengths[m] - (vs/60)*interval

                if lengths[m] < 0:

                    lengths[m] = 0
                    
                    states[m] = 0

    #########################################
    #########################################
                
        ###############
        #Plot length distribution every 5 minutes
        
        if i%300==0 and i != 0:
            
            fig, ax = plt.subplots(1,1, figsize=(4, 3.5))
            ax = plt.subplot(1,1,1)
            
            plt.hist(lengths, bins=np.arange(0,np.amax(lengths),np.amax(lengths)/15))
            
            plt.grid(False)
            ax.set_facecolor('white')
            ax.tick_params(axis='both', color = 'k', length = 4, width = 1.5)
            plt.xlim(0,)
            plt.ylim(0,)
            plt.title(str(int(i)) + " sec", fontsize = 15, color = 'k')
            plt.xlabel('Length (' + r"$\mu$" + 'm)', fontsize = 15, color = 'k')
            plt.ylabel('Frequency', fontsize = 15, color = 'k')
            plt.xticks(fontsize = 15, color = 'k')
            plt.yticks(fontsize = 15, color = 'k')
            plt.axhline(y = 0, color='k', linewidth=2.5)
            plt.axvline(x = 0, color='k', linewidth=2.5)
            #ax.text(0.1, 140, str(i*interval) + " sec", verticalalignment='center', horizontalalignment='left', fontsize = 15, color = "black")
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            #plt.savefig(savedir+str(int(i*interval))+'sec.jpg')
            plt.show()
            plt.tight_layout()
            
        ##############################

        lengthmt1.append(lengths[1])
        
        if sum(lengths) == 0:
            
            avglength = 0
            
        else:
        
            lengths = np.array(lengths)

            avglength = np.mean(lengths[np.nonzero(lengths)])

        avglengths_lst.append(avglength)
        
    ##############################

    fig, ax = plt.subplots(1,1, figsize=(10, 3.5))
    ax = plt.subplot(1,1,1)
    
    plt.plot(np.arange(0,tottime,interval), avglengths_lst, color = 'maroon', linewidth = 2)
    plt.plot(np.arange(0,tottime,interval), lengthmt1, color = 'grey', linewidth = 2, alpha = 0.65)
    
    plt.grid(False)
    ax.set_facecolor('white')
    ax.tick_params(axis='both', color = 'k', length = 4, width = 1.5)
    plt.xlim(0,)
    plt.ylim(0,)
    plt.xlabel('Time (sec)', fontsize = 15, color = 'k')
    plt.ylabel('Length (' + r"$\mu$" + 'm)', fontsize = 15, color = 'k')
    plt.xticks(fontsize = 15, color = 'k')
    plt.yticks(fontsize = 15, color = 'k')
    plt.axhline(y = 0, color='k', linewidth=2.5)
    plt.axvline(x = 0, color='k', linewidth=2.5)
    #ax.text(0.1, 140, str(i*interval) + " sec", verticalalignment='center', horizontalalignment='left', fontsize = 15, color = "black")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #plt.savefig(savedir+'Avglength over time.jpg')
    plt.show()
    plt.tight_layout()
    
    ##############################
    
    lifetimes = np.sort(np.array(lifetimes))
    
    cumfreq = []
    for i, d in enumerate(lifetimes):

        if i == 0:
            cumfreq.append(1.0/len(lifetimes))
        else:
            cumfreq.append((1.0/len(lifetimes)) + cumfreq[i-1])
    
    fig, ax = plt.subplots(1,1, figsize=(4, 3.5))
    ax = plt.subplot(1,1,1)
    
    plt.plot(lifetimes, cumfreq, color = 'k', linewidth=1.5)
    
    plt.plot(probrange, gammadis, linewidth = 2, color = 'r', linestyle = '--')
    
    plt.grid(False)
    ax.set_facecolor('white')
    ax.tick_params(axis='both', color = 'k', length = 4, width = 1.5)
    plt.xlim(0,900)
    plt.ylim(0,1)
    #plt.title('Lifetime OUT', fontsize = 15, color = 'k')
    plt.xlabel('Lifetime (sec)', fontsize = 15, color = 'k')
    plt.ylabel('Cumulative Probability', fontsize = 15, color = 'k')
    plt.xticks(fontsize = 15, color = 'k')
    plt.yticks(fontsize = 15, color = 'k')
    plt.axhline(y = 0, color='k', linewidth=2.5)
    plt.axvline(x = 0, color='k', linewidth=2.5)
    #ax.text(0.1, 140, str(i*interval) + " sec", verticalalignment='center', horizontalalignment='left', fontsize = 15, color = "black")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #plt.savefig(savedir+'Avglength over time.jpg')
    plt.show()
    plt.tight_layout()
    
    return(lifetimes)