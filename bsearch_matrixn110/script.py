#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 18-Dec-12


# This program controls the running of the Enzo simulator with an optimization routine of the input parameters
# This is a new version wihich creates multiple batches, each with multiple runs

from pylab import *
from subprocess import *
import sys
import os
import time

class Batch:
    def __init__(self,batchcounter):
        self.batchcounter=batchcounter
        self.done=False
        self.thisdir = Popen('pwd',shell=True, stdout = PIPE).communicate()[0].strip('\n')
        self.dd =self.thisdir+'/batchfiles/batch'+str(batchcounter)+'/'
        newdir=Popen('mkdir -p '+self.dd,shell=True)            # Create the new output directory                                               
        Popen.wait(newdir)
        copypoints=Popen('cp pointsfiles/points_'+str(batchcounter)+'.out '+self.dd+'points.out',shell=True)
        Popen.wait(copypoints)
        file=open(self.dd+'points.out')
        lines=file.readlines()
        file.close()
        self.NumPoints = len(lines)
        self.tasks = ['tran','enzo','enzo2','fom']
        for task in self.tasks:
            self.EditBatchFile(task)
            self.EditScriptFile(task)
        copyTS=Popen('cp transcript.csh '+self.dd,shell=True)
        Popen.wait(copyTS)
        return


    def EditBatchFile(self,task):
        filename =task+'batch'
        file=open(filename,'r')
        s=file.readlines()
        file.close()
        if task in ['enzo','enzo2']:
            scriptline = '    \"cd '+self.thisdir+' ;./'+task+'script.csh {};\"  \n\n'
            numline = '#PBS -l select='+str(8*self.NumPoints)+':ncpus=8:model=wes \n'
            runline = 'seq '+str(8*self.NumPoints)+' | parallel -j 1 -u --sshloginfile $PBS_NODEFILE\\\n'
        elif task == 'tran':
            numline = '#PBS -l select='+str(int(ceil(self.NumPoints/6.0)))+':ncpus=6:model=wes \n'
            runline = 'seq '+str(self.NumPoints)+' | parallel -j 6 -u --sshloginfile $PBS_NODEFILE\\\n'
            scriptline = '    \"cd '+self.thisdir+' ;./batchfiles/batch'+str(self.batchcounter)+'/ddfiles/run{}/transcript.csh {};\"  \n\n'
        elif task == 'fom':
            numline = '#PBS -l select='+str(int(ceil(self.NumPoints/6.0)))+':ncpus=6:model=wes \n'
            scriptline = 'cd '+self.thisdir+' \n'
            runline = '\n'
            s[10] = 'mpiexec -comm none -np '+str(self.NumPoints)+' fomscript.csh 0 \n'
        s[1] = numline
        s[7] = runline
        s[8] = scriptline
        file=open(filename,'w')
        for line in s:
            file.write(line)
        file.close()
        return

    def EditScriptFile(self,task):
        filename = task+'script.csh'
        file = open(filename,'r')
        s = file.readlines()
        file.close()
        if task in ['enzo','enzo2']:
            s[10] = '    cd ./batchfiles/batch'+str(self.batchcounter)+'/ddfiles/run$dir \n'
            s[30] = '    cd ./batchfiles/batch'+str(self.batchcounter)+'/ddfiles/run$dir \n'
            s[42] = '    cd ./batchfiles/batch'+str(self.batchcounter)+'/ddfiles/run$dir \n'
        elif task == 'fom':
            s[4] = 'cd ./batchfiles/batch'+str(self.batchcounter)+'/ddfiles/run$rank \n'
        else:
            s[3] = 'cd ./batchfiles/batch'+str(self.batchcounter)+'/ddfiles/run$dir \n'

        file=open(filename,'w')
        for line in s:
            file.write(line)
        file.close()
        return

    class Job:
            def __init__(self,batchcounter,counter):
                    self.counter=counter
                    self.batchcounter=batchcounter
                    self.EnzoBatchID = 0
                    self.done=False
                    thisdir = Popen('pwd',shell=True, stdout = PIPE).communicate()[0].strip('\n')
                    self.dd = thisdir+'/batchfiles/batch'+str(batchcounter)+'/ddfiles/run'+str(counter)+'/'
                    self.P = self.Parameters() # This holds all of the simulation parameters
                    self.P.ReadParams(counter,batchcounter) # Read the parameters from a file
                    EnzoWaitFile = "DD%04d/output_%04d"%(self.P.snapmax,self.P.snapmax)
                    self.WaitFiles=['ParticleVelocities.4','DD0005/output_0005','RunFinished','FomFinished','Done']

                    try:    # This checks to see if it is already done, for restarts
                            self.GetResult()
                            return
                    except IOError:
                            newdir=Popen('mkdir -p '+self.dd,shell=True) 		# Create the new output directory
                            Popen.wait(newdir)
                            #copySM1=Popen('cp smile.ini '+self.dd+'smilem.ini',shell=True) 	# Copy smilem.ini into new directory
                            #Popen.wait(copySM1)
                            #self.EditSmileIni(1)                                    # Edit smilem.ini
                            #copySS1=Popen('cp smilem.script '+self.dd,shell=True) 	# Copy smilem.script
                            #Popen.wait(copySS1)
                            #self.EditSmileScript(1)                                 # Edit smilem.script
                            #copySB1=Popen('cp smile.ini '+self.dd+'smileb.ini',shell=True) 	# Copy smileb.ini into new directory
                            #Popen.wait(copySB1)
                            #self.EditSmileIni(2)                                    # Edit smileb.ini
                            #copySS2=Popen('cp smileb.script '+self.dd,shell=True) 	# Copy smileb.script
                            #Popen.wait(copySS2)
                            #self.EditSmileScript(2)                                 # Edit smileb.script
                            copyTA1=Popen('cp triaxialm.cfg '+self.dd+'triaxialm.cfg',shell=True) 	# Copy triaxialm.cfg into new directory
                            Popen.wait(copyTA1)
                            self.EditTriaxialCfg(1)
                            copyTA2=Popen('cp triaxialb.cfg '+self.dd+'triaxialb.cfg',shell=True) 	# Copy triaxialm.cfg into new directory
                            Popen.wait(copyTA2)
                            self.EditTriaxialCfg(2)
                            copyTS=Popen('cp transcript.csh '+self.dd+'transcript.csh',shell=True) 	# Copy transcript.csh into new directory
                            Popen.wait(copyTS)
                            self.EditTranScript()
                            copyFI=Popen('cp fominput '+self.dd+'fominput',shell=True) 	# Copy fominput into new directory
                            Popen.wait(copyFI)
                            self.EditFomInput()
                            copyAM=Popen('cp AMRTest.enzo '+self.dd+'AMRTest.enzo',shell=True) 	# Copy AMRTest.enzo into new directory
                            Popen.wait(copyAM)
                            self.EditAMRTest()
                            self.waitfile = self.WaitFiles[0]
                    return


            class Parameters:
                def __init__(self):
                    self.ngasmain = 4000000 # Number of gas particles in the main cluster
                    self.nhalomain = 4000000 # Number of dm particles in the main cluster
                    self.ngasbullet = 1000000 # Number of gas particles in the bullet cluster
                    self.nhalobullet = 1000000 # Number of dm particles in the bullet cluster
                    self.rstart = 2800.0  # Starting separation
                    self.gtol=1.0E-8
                    self.fret=0.0
                    self.iter=0
                    self.snapmin=22
                    self.snapmax=30
                    self.NRotProcs = 12
                    self.timeoffirstsnap = 0.0
                    self.stepinGy = 0.9777
                    return

                def ReadParams(self,counter,batchcounter):
                    file=open('batchfiles/batch'+str(batchcounter)+'/points.out')
                    lines=file.readlines()
                    file.close()

                    for line in lines:
                            if int(line.strip().split()[0].split('=')[1]) == counter:
                                    self.C1=float(line.strip().split()[1].split('=')[1])
                                    self.C2=float(line.strip().split()[2].split('=')[1])
                                    self.M1=float(line.strip().split()[3].split('=')[1])
                                    self.M2=float(line.strip().split()[4].split('=')[1])
                                    self.RC1=float(line.strip().split()[5].split('=')[1])
                                    self.RC2=float(line.strip().split()[6].split('=')[1])
                                    self.Beta1=float(line.strip().split()[7].split('=')[1])
                                    self.Beta2=float(line.strip().split()[8].split('=')[1])
                                    self.Alpha1=float(line.strip().split()[9].split('=')[1])
                                    self.Alpha2=float(line.strip().split()[10].split('=')[1])
                                    self.Epsilon1=float(line.strip().split()[11].split('=')[1])
                                    self.Epsilon2=float(line.strip().split()[12].split('=')[1])
                                    self.Rs1=float(line.strip().split()[13].split('=')[1])
                                    self.Rs2=float(line.strip().split()[14].split('=')[1])
                                    self.Beta21=float(line.strip().split()[15].split('=')[1])
                                    self.Beta22=float(line.strip().split()[16].split('=')[1])
                                    self.RC21=float(line.strip().split()[17].split('=')[1])
                                    self.RC22=float(line.strip().split()[18].split('=')[1])
                                    self.N21=float(line.strip().split()[19].split('=')[1])
                                    self.N22=float(line.strip().split()[20].split('=')[1])
                                    self.GF1=float(line.strip().split()[21].split('=')[1])
                                    self.GF2=float(line.strip().split()[22].split('=')[1])
                                    self.Z=float(line.strip().split()[23].split('=')[1]) 
                                    self.Mag=float(line.strip().split()[24].split('=')[1]) 
                                    self.FieldMode=line.strip().split()[25].split('=')[1]
                                    self.Tbg=float(line.strip().split()[26].split('=')[1]) 
                                    self.LogRhoMin=float(line.strip().split()[27].split('=')[1]) 
                                    self.IP=float(line.strip().split()[28].split('=')[1]) 
                                    self.IPPhi=float(line.strip().split()[29].split('=')[1]) 
                                    self.IPTheta=float(line.strip().split()[30].split('=')[1]) 
                                    self.IPPsi=float(line.strip().split()[31].split('=')[1]) 
                                    self.Vinc=float(line.strip().split()[32].split('=')[1]) 
                                    self.P1=float(line.strip().split()[33].split('=')[1])# These control the main shape 
                                    self.Q1=float(line.strip().split()[34].split('=')[1]) 
                                    self.P2=float(line.strip().split()[35].split('=')[1])# These control the bullet shape 
                                    self.Q2=float(line.strip().split()[36].split('=')[1]) 
                                    self.phi1=float(line.strip().split()[37].split('=')[1])# These control the main orientation
                                    self.theta1=float(line.strip().split()[38].split('=')[1])
                                    self.psi1=float(line.strip().split()[39].split('=')[1]) 
                                    self.phi2=float(line.strip().split()[40].split('=')[1]) # These control the bullet orientation
                                    self.theta2=float(line.strip().split()[41].split('=')[1]) 
                                    self.psi2=float(line.strip().split()[42].split('=')[1])  
                                    self.Visc=float(line.strip().split()[43].split('=')[1])  
                                    self.TFudge=float(line.strip().split()[44].split('=')[1])  
                                    break	
                            else:
                                    continue		
                    self.Phi=0.0 # These control the viewing angle
                    self.Theta=1.11
                    self.ThetaMin = 0.0
                    self.ThetaMax = 3.14
                    self.Psi=3.03
                    self.PsiMin = 0.0
                    self.PsiMax = 6.28
                    self.Vx=4336.0
                    self.Vy=728.0
                    self.Vz=-1315.0
                    return


            def WaitFileExists(self):
                    # This subroutine checks to see if the waitfile exists and is complete
                    NumTries = 0
                    LastFileSize = 0
                    #print "Looking for WaitFile ",self.dd+self.waitfile
                    #sys.stdout.flush()
                    while NumTries < 4:
                            try:
                                    NumTries = NumTries + 1
                                    FileSize = os.path.getsize(self.dd+self.waitfile)
                                    #print "WaitFile = %s, FileSize = %d, LastFileSize = %d\n"%(self.waitfile,FileSize,LastFileSize)
                                    #sys.stdout.flush()
                                    if FileSize > 0 and FileSize == LastFileSize :
                                        #print "Returning True"
                                        #sys.stdout.flush()
                                        return True
                                    else:
                                        LastFileSize = FileSize
                                        time.sleep(1.0)
                                        continue
                            except OSError:
                                time.sleep(1.0)
                    #print "Returning False"
                    #print "WaitFile = %s, FileSize = %d, LastFileSize = %d\n"%(self.dd+self.waitfile,FileSize,LastFileSize)
                    #sys.stdout.flush()
                    #sys.exit()
                    return False			


            def GetResult(self):
                # This finds the result of the FindFom run.
                nfom = 'newfom_massx1.out' # This identifies the FOM file
                file = open(self.dd+nfom,'r')
                line=file.readline()
                file.close()
                self.fom=float(line.strip().split()[0].split('FOM=')[1].split(',')[0])
                self.bestsnap=int(line.strip().split()[5].split(',')[0])
                return


            def EditTranScript(self):
                    filename = self.dd+'transcript.csh' 
                    file=open(filename,'r')
                    s=file.readlines()
                    file.close()
                    s[15]='/home1/clage/Research/bullet/code/CombineGalaxies_NewerIP/CombineGalaxies bullet.dat '+str(self.P.phi2)+' '+str(self.P.theta2)+' '+str(self.P.psi2)+' main.dat '+str(self.P.phi1)+' '+str(self.P.theta1)+' '+str(self.P.psi1)+' '+str(self.P.IP)+' '+str(self.P.IPPhi)+' '+str(self.P.IPTheta)+' '+str(self.P.IPPsi)+' '+str(self.P.rstart)+' '+str(self.P.Vinc)+' collision.dat'+' >&collision.out \n'# New version of CombineGalaxies 
                    s[19]='python /home1/clage/Research/bullet/code/pysubs/make_cool_tfudge.py '+str(self.P.Z)+' '+str(self.P.TFudge)+' >&cool.out \n'
                    s[21]='python /home1/clage/Research/bullet/code/pysubs/translate_gadget_to_enzo_batch_15Jun13.py '+str(self.P.Mag)+' '+str(self.P.FieldMode)+' '+str(self.P.Tbg)+' '+str(self.P.LogRhoMin)+' >&tran.out \n'
                    file=open(filename,'w')
                    for line in s:
                            file.write(line)
                    file.close()
                    return

            def EditTriaxialCfg(self,Cluster): 
                # Edits the triaxial.cfg file
                # 1 = Main, 2 = Bullet
                if Cluster == 1:
                    filename =self.dd+'triaxialm.cfg' 
                    file=open(filename,'r')
                    s=file.readlines()
                    file.close()
                    s[1]='P = '+str(self.P.P1)+'\n'
                    s[2]='Q = '+str(self.P.Q1)+'\n'
                    s[4]='M200 = '+str(self.P.M1)+'\n'
                    s[5]='C = '+str(self.P.C1)+'\n'
                    s[9]='NDM = '+str(self.P.nhalomain)+'\n'
                    s[10]='NGas = '+str(self.P.ngasmain)+'\n'
                    s[11]='GF = '+str(self.P.GF1)+'\n'
                    s[13]='GasAlpha = '+str(self.P.Alpha1)+'\n'
                    s[14]='GasBeta = '+str(self.P.Beta1)+'\n'
                    s[15]='GasRcool = '+str(self.P.RC1)+'\n'
                    s[16]='GasEpsilon = '+str(self.P.Epsilon1)+'\n'
                    s[17]='GasRs = '+str(self.P.Rs1)+'\n'
                    s[18]='GasBeta2 = '+str(self.P.Beta21)+'\n'
                    s[19]='GasRcool2 = '+str(self.P.RC21)+'\n'
                    s[20]='GasN2 = '+str(self.P.N21)+'\n'
                    s[22] = 'potentialfile 	= ../../../../potentialm.dat\n'
                    s[23] = 'nbodyfile	= ../../../../nbodym.dat\n'
                elif Cluster == 2:
                    filename =self.dd+'triaxialb.cfg' 
                    file=open(filename,'r')
                    s=file.readlines()
                    file.close()
                    s[1]='P = '+str(self.P.P2)+'\n'
                    s[2]='Q = '+str(self.P.Q2)+'\n'
                    s[4]='M200 = '+str(self.P.M2)+'\n'
                    s[5]='C = '+str(self.P.C2)+'\n'
                    s[9]='NDM = '+str(self.P.nhalobullet)+'\n'
                    s[10]='NGas = '+str(self.P.ngasbullet)+'\n'
                    s[11]='GF = '+str(self.P.GF2)+'\n'
                    s[13]='GasAlpha = '+str(self.P.Alpha2)+'\n'
                    s[14]='GasBeta = '+str(self.P.Beta2)+'\n'
                    s[15]='GasRcool = '+str(self.P.RC2)+'\n'
                    s[16]='GasEpsilon = '+str(self.P.Epsilon2)+'\n'
                    s[17]='GasRs = '+str(self.P.Rs2)+'\n'
                    s[18]='GasBeta2 = '+str(self.P.Beta22)+'\n'
                    s[19]='GasRcool2 = '+str(self.P.RC22)+'\n'
                    s[20]='GasN2 = '+str(self.P.N22)+'\n'
                    s[22] = 'potentialfile 	= ../../../../potentialb.dat\n'
                    s[23] = 'nbodyfile	= ../../../../nbodyb.dat\n'
                file=open(filename,'w')
                for line in s:
                        file.write(line)
                file.close()
                return

            def EditSmileIni(self,Cluster): 
                # Edits the smile.ini file
                # 1 = Main, 2 = Bullet
                if Cluster == 1:
                    filename =self.dd+'smilem.ini' 
                    file=open(filename,'r')
                    s=file.readlines()
                    file.close()
                    s[5]='P = '+str(self.P.P1)+'\n'
                    s[6]='Q = '+str(self.P.Q1)+'\n'
                    s[8]='Rc = '+str(self.P.C1)+'\n'
                elif Cluster == 2:
                    filename =self.dd+'smileb.ini' 
                    file=open(filename,'r')
                    s=file.readlines()
                    file.close()
                    s[5]='P = '+str(self.P.P2)+'\n'
                    s[6]='Q = '+str(self.P.Q2)+'\n'
                    s[8]='Rc = '+str(self.P.C2)+'\n'
                file=open(filename,'w')
                for line in s:
                        file.write(line)
                file.close()
                return

            def EditSmileScript(self,Cluster): 
                # Edits the smile.script file
                # 1 = Main, 2 = Bullet
                if Cluster == 1:
                    filename =self.dd+'smilem.script' 
                    file=open(filename,'r')
                    s=file.readlines()
                    file.close()
                    s[4]='ExportNbody('+str(self.P.nhalomain)+',"nbodym.dat")\n'
                elif Cluster == 2:
                    filename =self.dd+'smileb.script' 
                    file=open(filename,'r')
                    s=file.readlines()
                    file.close()
                    s[4]='ExportNbody('+str(self.P.nhalobullet)+',"nbodyb.dat")\n'
                file=open(filename,'w')
                for line in s:
                        file.write(line)
                file.close()
                return

            def EditAMRTest(self): 
                # Edits the AMRTest.enzo file
                filename =self.dd+'AMRTest.enzo' 
                file=open(filename,'r')
                s=file.readlines()
                file.close()
                DensityUnits = 6.76976638e-22
                smallrho = 10**self.P.LogRhoMin / 3.0 * DensityUnits
                s[94]='SmallRho = '+str(smallrho)+'\n'
                s[100]='ViscosityCoefficient = '+str(self.P.Visc)+'\n'
                file=open(filename,'w')
                for line in s:
                        file.write(line)
                file.close()
                return

            def EditFomInput(self): 
                # Edits the fominput file
                filename =self.dd+'fominput' 
                file=open(filename,'r')
                s=file.readlines()
                file.close()
                s[0] = 'Z '+str(self.P.Z)+' \n'
                s[1] = 'snapmin '+str(self.P.snapmin)+' \n'
                s[2] = 'snapmax '+str(self.P.snapmax)+' \n'
                s[3] = 'Theta '+str(self.P.Theta)+' \n'
                s[4] = 'Psi '+str(self.P.Psi)+' \n'
                s[5] = 'Vx '+str(self.P.Vx)+' \n'
                s[6] = 'Vy '+str(self.P.Vy)+' \n'
                s[7] = 'Vz '+str(self.P.Vz)+' \n'
                s[8] = 'TFudge '+str(self.P.TFudge)+' \n'
                file=open(filename,'w')
                for line in s:
                        file.write(line)
                file.close()
                return

            def CleanUp(self): 
                # Deletes files that are no longer needed
                print "Cleaning up files in directory %s\n"%self.dd

                delf = Popen('rm -f '+self.dd+'Grid*',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f '+self.dd+'Par*',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f '+self.dd+'Tot*',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f '+self.dd+'Int*',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f '+self.dd+'*.dat',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f  '+self.dd+'*.nb',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f  '+self.dd+'fom.out',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f  '+self.dd+'perf*',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f  '+self.dd+'enzo.out',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f  '+self.dd+'core*',shell=True)
                Popen.wait(delf)
                delf = Popen('rm -f  '+self.dd+'enzo.out',shell=True)
                Popen.wait(delf)
                keep = []#[max(self.bestsnap - 1, self.P.snapmin), self.bestsnap, min(self.bestsnap+1, self.P.snapmax)]
                for i in range(132):
                    if i in keep:
                        continue
                    else:
                        delf = Popen('rm -rf '+self.dd+'DD%.4d'%i,shell=True)
                        Popen.wait(delf)
                return


    def RunBatch(self):
        Jobs = range(1,self.NumPoints+1)
        PopJobs = list(Jobs)
        PopJobs.reverse()
        FinishedJobs=[]
        ActiveJobs=[]
        MaxNumberinParallel = self.NumPoints
        while len(ActiveJobs) < MaxNumberinParallel:
            newjob = PopJobs.pop()
            ActiveJobs.append(self.Job(self.batchcounter,newjob))
            print 'Starting job # %d\n'%newjob
        print 'All jobs started'

        Commands=['qsub -q devel -N tran'+str(self.batchcounter)+' tranbatch','qsub -q devel -N enzo_'+str(self.batchcounter)+' enzobatch','qsub -q devel -N enzo2_'+str(self.batchcounter)+' enzo2batch', 'qsub -q devel -N fom'+str(self.batchcounter)+' fombatch']

        for i, command in enumerate(Commands):
            while True: # Keep trying until job can be launched successfully
                try:
                    BatchJobID = int(Popen(command, shell=True, stdout=PIPE).communicate()[0].split('.')[0])
                    break
                except:
                    time.sleep(15.0)

            BatchJobStateCheck = 'qstat -fx %d |grep job_state'%BatchJobID
            JobsStillRunning = True
            JobsCompletedSuccessfully = False
            for Job in ActiveJobs:
                Job.waitfile = Job.WaitFiles[i]

            while JobsStillRunning:
                time.sleep(2.0)
                BatchJobState = Popen(BatchJobStateCheck, shell=True, stdout=PIPE).communicate()[0].split()[2]
                for Job in ActiveJobs:
                    if Job.WaitFileExists():
                        JobsStillRunning = False
                    else:
                        JobsStillRunning = True
                        break
                if BatchJobState == 'F' and not JobsStillRunning:
                    JobsCompletedSuccessfully = True
                elif BatchJobState == 'F' and JobsStillRunning:
                    JobsStillRunning = False
                    JobsCompletedSuccessfully = False

        #if JobsCompletedSuccessfully:
        for Job in ActiveJobs:
            Job.CleanUp()
        return


def roundup(i):
    # Rounds up the number of points in the points.out file until it is a multiple of 12
    k = i / 12
    if i%12 >0:
        k = k+1
    return k * 12


def build_delta_points(batchcounter, names, bestbatch, bestcounter):
    # Builds a points files with random deltas
    deltas = [0.95, 1.05]
    outfile = open('fill_space.cfg','w')
    for name in names:
        outfile.write(name+' 1.0 '+str(deltas[0])+' '+str(deltas[1])+'\n')
    outfile.close()
    fillbatchcmd = 'python ~/Research/bullet/code/pysubs/fill_space/fill_space4.py batchfiles/batch'+str(bestbatch)+'/points_massx1.out '+str(bestcounter)+' 32 5 fill_space.cfg pointsfiles/points_'+str(batchcounter)+'.out'
    print fillbatchcmd

    fillbatch=Popen(fillbatchcmd,shell=True)
    # Now wait for output file to be done
    NumTries = 0
    LastFileSize = 0
    while NumTries < 50:
        try:
            NumTries = NumTries + 1
            FileSize = os.path.getsize('pointsfiles/points_'+str(batchcounter)+'.out')
            print FileSize
            if FileSize > 0 and FileSize == LastFileSize :
                return
            else:
                LastFileSize = FileSize
                time.sleep(10.0)
                continue
        except OSError:
            print "File not there"
            time.sleep(10.0)
    return

def make_combined_file(batchcounter):

    outfile=open('batchfiles/batch'+str(batchcounter)+'/points_massx1.out','w')
    fom=list()
    infile=open('batchfiles/batch'+str(batchcounter)+'/points.out')
    lines=infile.readlines()
    infile.close()
    numpoints = len(lines)

    for run in range(1,numpoints+1):

            for pointsline in lines:
                    pointsrun = int(pointsline.strip().split()[0].split('=')[1])
                    if pointsrun == run:
                            break
            fomfile = 'newfom_massx1.out'

            infile = 'batchfiles/batch'+str(batchcounter)+'/ddfiles/run'+str(run)+'/'+fomfile
            try:
                    data = open(infile,'r')
                    line = data.readline()
                    fom = float(line.strip().split()[0].strip('FOM=').strip(','))
                    xfom = float(line.strip().split()[1].strip('XFOM=').strip(','))
            except:
                    fom = -1.0
                    xfom = -1.0
            print numpoints, run, fom, xfom
            pointsline = pointsline.strip('\n')
            pointsline = pointsline+' Chi2='+str(fom)+' XChi2='+str(xfom)+'\n'
            pointsline = pointsline
            outfile.write(pointsline)

    outfile.close()
    return

def read_combined_files(batchcounter, fullnames):
    bestchi2 = 1000.0
    points=[]
    values=[]
    for i in range(batchcounter):
        infilename = 'batchfiles/batch'+str(i)+'/points_massx1.out'
        file = open(infilename,'r')
        lines = file.readlines()
        file.close()

        NumEval  = len(lines)
        NumParam = len(fullnames)
        print "NumEval = %d, NumParam = %d\n"%(NumEval,NumParam)
        for line in lines:
            point = []
            entries = line.strip().split()
            for entry in entries:
                name = entry.split('=')[0]
                if name == "Counter":
                    counter = int(entry.split('=')[1])
                if name in fullnames:
                    if name == 'FieldMode':
                        point.append(entry.split('=')[1])
                    else:
                        point.append(float(entry.split('=')[1]))

                if name =='Chi2':
                    chi2 = float(entry.split('=')[1])
            points.append(point)
            values.append(chi2)
            if chi2 >0 and chi2 < bestchi2:
                bestbatch = i
                bestchi2 = chi2
                bestcounter = counter
                bestindex = len(values) - 1
        print "bestbatch = %d, bestcounter = %d"%(bestbatch,bestcounter)
    return [points,values,bestindex,bestbatch,bestcounter]


def write_log_file(batchcounter,fullnames,scale,bestchi2):
    logfilename = 'script.log'
    logfile = open(logfilename, 'a')
    logfile.write('Counter=%d '%batchcounter)
    for i,name in enumerate(fullnames):
        logfile.write(name+'='+str(scale[i])+' ')
    logfile.write('\n')
    logfile.write('BestChi2 = %.4f\n'%bestchi2)
    logfile.write('\n')
    logfile.close()
    return
                      

#****************MAIN PROGRAM*****************
# multi_lowresn81 run 96
# Counter=96 C1=1.5857 C2=7.1550 M1=260863.8786 M2=29193.2792 RC1=68.8907 RC2=46.8176 Beta1=0.4176 Beta2=1.0067 Alpha1=0.7513 Alpha2=1.1442 Epsilon1=0.0000 Epsilon2=0.0000 Rs1=19.9867 Rs2=18.0567 Beta21=0.7346 Beta22=0.4436 RC21=644.1680 RC22=507.2446 N21=0.0000 N22=0.0000 GF1=0.1750 GF2=0.1618 Z=0.7853 Mag=7.2932 FieldMode=RandomRadial  Tbg=4.0120 LogRhoMin=-8.0000 IP=248.6520 IPTheta=94.7669 Vinc=0.9874 P1=0.3411 Q1=0.6800 P2=0.6211 Q2=0.6837 Phi1=184.2472 Theta1=39.2320 Psi1=199.1037 Phi2=165.6049 Theta2=101.4386 Psi2=64.9385 Visc=0.1375 Chi2=1.628145

fullnames=['C1','C2', 'M1', 'M2', 'RC1','RC2','Beta1', 'Beta2', 'Alpha1', 'Alpha2', 'Epsilon1', 'Epsilon2', 'Rs1', 'Rs2', 'Beta21', 'Beta22', 'RC21', 'RC22', 'N21', 'N22', 'GF1','GF2','Z','Mag',  'FieldMode', 'Tbg', 'LogRhoMin', 'IP', 'IPPhi', 'IPTheta', 'IPPsi', 'Vinc', 'P1', 'Q1', 'P2', 'Q2', 'Phi1', 'Theta1', 'Psi1', 'Phi2', 'Theta2', 'Psi2', 'Visc', 'TFudge']

names=['RC1','RC2','Beta1', 'Beta2', 'Alpha1', 'Alpha2', 'Rs1', 'Rs2', 'Beta21', 'Beta22', 'RC21', 'RC22', 'GF1','GF2', 'Mag', 'Vinc', 'Visc']

#minvalues=[1.0, 3.0, 100000.0, 10000.0, 10.00, 10.00, 0.20, 0.20, 0.05, 0.05, 0.00, 0.00, 1.0, 1.0, 0.20, 0.20, 50.0, 50.0, 0.00, 0.00, 0.10, 0.10, 0.10, 1.00, 'RandomRadial', 1.0, -9.0, 0.0, 0.0, 0.0, 0.0, 0.80, 0.30, 0.30, 0.30, 0.30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.3]

#maxvalues=[4.0, 12.0, 350000.0, 35000.0, 500.00, 500.00, 1.0, 1.5, 1.5, 1.5, 0.00, 0.00, 100.0, 100.0, 2.0, 2.0, 2000.0, 2000.0, 0.00, 0.00, 0.28, 0.28, 0.99, 100.0, 'RandomRadial', 8.0, -5.0, 300.0, 360.0, 180.0, 360.0, 1.20, 0.99, 0.99, 0.99, 0.99, 360.0, 180.0, 360.0, 360.0, 180.0, 360.0, 0.50, 1.0]

Iterations = 10 # of iterations of iterations of gradient finding.
for batchcounter in range(2, Iterations+1):
    [points,values,bestindex,bestbatch,bestcounter] = read_combined_files(batchcounter, fullnames)
    print "bestbatch = %d, bestcounter = %d"%(bestbatch,bestcounter)
    scale = points[bestindex] # Best value so far
    bestchi2 = values[bestindex]
    write_log_file(batchcounter,fullnames,scale,bestchi2)
    build_delta_points(batchcounter, names, bestbatch, bestcounter)
    newbatch = Batch(batchcounter)
    newbatch.RunBatch() # Run the delta points
    make_combined_file(batchcounter)
    

#************END MAIN PROGRAM*************************


