import numpy as np
from Vepp_lib import freq_lib as fl
class Twiss_data():
    '''
    BF-бета из фазы
    BA-бета из амплитуды
    exp_data=[f1,f2,f3,freq]
    model_data=importData
    фаза 4 по счету
    '''
    def __init__(self,type:str,bpm_numbers:np.array,x_coordinats1:np.array,x_coordinats2:np.array,x_coordinats3:np.array,model_phase:np.array,model_beta:np.array,model_alf:np.array,len:int):

        self.type=type
        self.x_coordinats1=x_coordinats1
        self.x_coordinats2 = x_coordinats2
        self.x_coordinats3 = x_coordinats3
        self.model_alf=model_alf
        self.model_phase = model_phase
        self.model_beta=model_beta
        self.len=len
        self.bpm_numbers=bpm_numbers

    def beta_from_ampl(self,number,all_data_massive,len):
        self.class_massive=np.empty(0)
        self.ampl_massive=np.empty(0)
        for i in range(54):
            self.class_massive=np.append(self.class_massive,fl.Calculate_Freq(all_data_massive[i][:len]))
            self.ampl_massive=np.append(self.ampl_massive,self.class_massive[i].calculate_ampl_and_phase(len)[0])

        self.action=np.mean(self.ampl_massive**2/(2*self.model_beta[:54]))
        self.beta=self.ampl_massive[number]**2/(2*self.action)
        return self.beta

    def exp_phase(self,exp_ampl_and_phase,model_delt_bpm_phase,parametr):
        self.exp_delt_phase = np.empty(0)
        for i in np.arange(parametr).tolist():
            self.delt_exp_phase_massive = np.array(
                [np.abs(exp_ampl_and_phase[i + 1][1]) - np.abs(exp_ampl_and_phase[i][1]),
                 np.abs(exp_ampl_and_phase[i + 1][1]) + np.abs(exp_ampl_and_phase[i][1]),
                 -np.abs(exp_ampl_and_phase[i + 1][1]) + np.abs(exp_ampl_and_phase[i][1]),
                 -np.abs(exp_ampl_and_phase[i + 1][1]) - np.abs(exp_ampl_and_phase[i][1])])

            self.compare_values = np.abs(
                np.array([1 / np.tan(self.delt_exp_phase_massive[0]) - 1 / np.tan(model_delt_bpm_phase[i]),
                          1 / np.tan(self.delt_exp_phase_massive[1]) - 1 / np.tan(model_delt_bpm_phase[i]),
                          1 / np.tan(self.delt_exp_phase_massive[2]) - 1 / np.tan(model_delt_bpm_phase[i]),
                          1 / np.tan(self.delt_exp_phase_massive[3]) - 1 / np.tan(model_delt_bpm_phase[i])]))
            self.index = np.argmin(self.compare_values)
            self.exp_delt_phase = np.append(self.exp_delt_phase, self.delt_exp_phase_massive[self.index])
        return self.exp_delt_phase

    def beta_from_phase_3(self,x_coordinats1,x_coordinats2,x_coordinats3,bpm_numbers,type):

        self.x_coordinats1=x_coordinats1
        self.x_coordinats2 = x_coordinats2
        self.x_coordinats3 = x_coordinats3
        self.bpm_numbers=bpm_numbers

        '''
        пока что для BF используются 3 соседних пикапа
        '''

        self.parametrs=np.array([fl.Calculate_Freq(self.x_coordinats1),fl.Calculate_Freq(self.x_coordinats2),fl.Calculate_Freq(self.x_coordinats3)])
        self.freq=self.parametrs[1].calculate_freq
        if self.freq < 0.5:
            self.freq=1-self.freq
        #self.freq = 0.5368830987373664
        self.exp_ampl_and_phase=np.array([self.parametrs[0].calculate_ampl_and_phase(self.len),self.parametrs[1].calculate_ampl_and_phase(self.len),self.parametrs[2].calculate_ampl_and_phase(self.len)])

        if (self.bpm_numbers[0]<self.bpm_numbers[1]<self.bpm_numbers[2]):

            self.model_delt_bpm_phase = np.array(
                [self.model_phase[self.bpm_numbers[1]] - self.model_phase[self.bpm_numbers[0]],
                 self.model_phase[self.bpm_numbers[2]] - self.model_phase[self.bpm_numbers[1]]])

            self.exp_delt_phase = self.exp_phase(self.exp_ampl_and_phase, self.model_delt_bpm_phase, 2)

        if (self.bpm_numbers[0]>self.bpm_numbers[1]>self.bpm_numbers[2]):

            self.model_delt_bpm_phase = np.array(
                [self.model_phase[self.bpm_numbers[1]] - self.model_phase[self.bpm_numbers[2]],
                 self.model_phase[self.bpm_numbers[0]] - self.model_phase[self.bpm_numbers[1]]])

            self.exp_delt_phase = self.exp_phase(np.flipud(self.exp_ampl_and_phase), self.model_delt_bpm_phase, 2)



        if (self.bpm_numbers[0] > self.bpm_numbers[1] < self.bpm_numbers[2]):

            self.model_delt_bpm_phase = np.array(
                [ self.model_phase[self.bpm_numbers[0]]-self.model_phase[self.bpm_numbers[1]],
                 self.model_phase[self.bpm_numbers[2]] - self.model_phase[self.bpm_numbers[1]]])


            self.exp_delt_phase0 =2*np.pi*self.freq-self.exp_phase(
                np.array([self.exp_ampl_and_phase[0], self.exp_ampl_and_phase[1], self.exp_ampl_and_phase[0]]),
                self.model_delt_bpm_phase, 1)[0]

            self.exp_delt_phase1 = self.exp_phase(
                np.array([self.exp_ampl_and_phase[2], self.exp_ampl_and_phase[1], self.exp_ampl_and_phase[2]]),
                self.model_delt_bpm_phase, 2)[1]

            self.exp_delt_phase=np.array([self.exp_delt_phase0,self.exp_delt_phase1])
            self.model_delt_bpm_phase=np.array([self.model_phase[54]-self.model_delt_bpm_phase[0],self.model_delt_bpm_phase[1]])
        if (self.bpm_numbers[0] < self.bpm_numbers[1] > self.bpm_numbers[2]):

            self.model_delt_bpm_phase = np.array(
                [ self.model_phase[self.bpm_numbers[1]] - self.model_phase[self.bpm_numbers[0]],
                  self.model_phase[self.bpm_numbers[1]] - self.model_phase[self.bpm_numbers[2]]])

            self.exp_delt_phase0 = self.exp_phase(
                np.array([ self.exp_ampl_and_phase[1], self.exp_ampl_and_phase[0],
                          self.exp_ampl_and_phase[1]]),
                self.model_delt_bpm_phase, 1)[0]

            self.exp_delt_phase1 = 2*np.pi*self.freq-self.exp_phase(
                np.array([self.exp_ampl_and_phase[1], self.exp_ampl_and_phase[2], self.exp_ampl_and_phase[1]]),
                self.model_delt_bpm_phase, 2)[1]
            self.exp_delt_phase = np.array([self.exp_delt_phase0, self.exp_delt_phase1])

            self.model_delt_bpm_phase=np.array([self.model_delt_bpm_phase[0],self.model_phase[54]-self.model_delt_bpm_phase[1]])
        if type==1:
            self.beta = self.model_beta[self.bpm_numbers[1]] * (
                    (1 / np.tan(self.exp_delt_phase[0]) + 1 / np.tan(self.exp_delt_phase[1])) / (
                    1 / np.tan(self.model_delt_bpm_phase[0]) + 1 / np.tan(self.model_delt_bpm_phase[1])))
        if type==0:
            self.beta = self.model_beta[self.bpm_numbers[0]] * (
                    (1 / np.tan(self.exp_delt_phase[0]) - 1 / np.tan(self.exp_delt_phase[1]+self.exp_delt_phase[0])) / (
                    1 / np.tan(self.model_delt_bpm_phase[0]) - 1 / np.tan(self.model_delt_bpm_phase[1]+self.model_delt_bpm_phase[0])))
        if type==2:
            self.beta = self.model_beta[self.bpm_numbers[2]] * (
                    (1 / np.tan(self.exp_delt_phase[1]) - 1 / np.tan(self.exp_delt_phase[1]+self.exp_delt_phase[0])) / (
                    1 / np.tan(self.model_delt_bpm_phase[1]) - 1 / np.tan(self.model_delt_bpm_phase[1]+self.model_delt_bpm_phase[0])))
        if type==10:
            self.beta=self.model_alf[self.bpm_numbers[0]]*(
                    (1 / np.tan(self.exp_delt_phase[0]) - 1 / np.tan(self.exp_delt_phase[1]+self.exp_delt_phase[0])) / (
                    1 / np.tan(self.model_delt_bpm_phase[0]) - 1 / np.tan(self.model_delt_bpm_phase[1]+self.model_delt_bpm_phase[0])))+(1 / np.tan(self.exp_delt_phase[0])*1 / np.tan(self.model_delt_bpm_phase[0]+self.model_delt_bpm_phase[1])-1 / np.tan(self.model_delt_bpm_phase[0])*1 / np.tan(self.exp_delt_phase[1]+self.exp_delt_phase[0]))/(1 / np.tan(self.model_delt_bpm_phase[0]) - 1 / np.tan(self.model_delt_bpm_phase[1]+self.model_delt_bpm_phase[0]))
        if type==11:
            self.beta=self.model_alf[self.bpm_numbers[1]]*(
                    (1 / np.tan(self.exp_delt_phase[0]) + 1 / np.tan(self.exp_delt_phase[1])) / (
                    1 / np.tan(self.model_delt_bpm_phase[0]) + 1 / np.tan(self.model_delt_bpm_phase[1])))+(1 / np.tan(self.exp_delt_phase[0])*1 / np.tan(self.model_delt_bpm_phase[1])-1 / np.tan(self.model_delt_bpm_phase[0])*1 / np.tan(self.exp_delt_phase[1]))/(1 / np.tan(self.model_delt_bpm_phase[0]) + 1 / np.tan(self.model_delt_bpm_phase[1]))

        if type==12:
            self.beta=self.model_alf[self.bpm_numbers[2]]*(
                    (1 / np.tan(self.exp_delt_phase[0]+self.exp_delt_phase[1]) - 1 / np.tan(self.exp_delt_phase[1])) / (
                    1 / np.tan(self.model_delt_bpm_phase[0]+self.model_delt_bpm_phase[1]) - 1 / np.tan(self.model_delt_bpm_phase[1])))+(1 / np.tan(self.exp_delt_phase[1])*1 / np.tan(self.model_delt_bpm_phase[0]+self.model_delt_bpm_phase[1])-1 / np.tan(self.model_delt_bpm_phase[1])*1 / np.tan(self.exp_delt_phase[1]+self.exp_delt_phase[0]))/(1 / np.tan(self.model_delt_bpm_phase[0]+self.model_delt_bpm_phase[1]) - 1 / np.tan(self.model_delt_bpm_phase[1]))
        return self.beta

    def dphase(self):
        return self.exp_delt_phase

    def beta_from_phase_N(self,N,number,data):

        self.coordinats_massive=np.empty([0])
        self.beta_massive1=np.empty([0])
        self.beta_massive0 = np.empty([0])
        self.beta_massive2 = np.empty([0])
        self.alpha_massive0=np.empty([0])
        self.alpha_massive1 = np.empty([0])
        self.alpha_massive2 = np.empty([0])
        self.index_BPM=np.empty([0])
        for i in np.arange(2*N+1).tolist():
            self.flag=number-N+i
            if self.flag<0:
                self.flag=self.flag+54
            if self.flag>53:
                self.flag = self.flag - 54
            self.coordinats_massive=np.append(self.coordinats_massive,Trajectory(init_cond=np.array([[1,0]]),sext_location=3,strenght=0.0,number=self.flag,power=4024,data=data).calculate_trajectory()[0])
            self.index_BPM=np.append(self.index_BPM,self.flag)
        self.coordinats_massive=np.reshape(self.coordinats_massive,(2*N+1,4024))
        '''
        0-крайний левый
        1-средний
        2-крайний правый

        '''
        for i in np.arange(1,N+1).tolist():
            for j in np.arange(1,N+1).tolist():
                self.beta_massive1=np.append(self.beta_massive1,self.beta_from_phase_3(self.coordinats_massive[N-i],self.coordinats_massive[N],self.coordinats_massive[N+j],np.array([int(self.index_BPM[N-i]),int(self.index_BPM[N]),int(self.index_BPM[N+j])]),1))
                self.alpha_massive1=np.append(self.alpha_massive1,self.beta_from_phase_3(self.coordinats_massive[N-i],self.coordinats_massive[N],self.coordinats_massive[N+j],np.array([int(self.index_BPM[N-i]),int(self.index_BPM[N]),int(self.index_BPM[N+j])]),11))
        for i in np.arange(N+1,2*N).tolist():
            for j in np.arange(i+1,2*N+1).tolist():
                self.beta_massive0=np.append(self.beta_massive0,self.beta_from_phase_3(self.coordinats_massive[N],self.coordinats_massive[i],self.coordinats_massive[j],np.array([int(self.index_BPM[N]),int(self.index_BPM[i]),int(self.index_BPM[j])]),0))
                self.alpha_massive0=np.append(self.alpha_massive0,self.beta_from_phase_3(self.coordinats_massive[N],self.coordinats_massive[i],self.coordinats_massive[j],np.array([int(self.index_BPM[N]),int(self.index_BPM[i]),int(self.index_BPM[j])]),10))
        for i in np.arange(N-1,0,-1).tolist():
            for j in np.arange(i-1,-1,-1).tolist():
                self.beta_massive2=np.append(self.beta_massive2,self.beta_from_phase_3(self.coordinats_massive[j],self.coordinats_massive[i],self.coordinats_massive[N],np.array([int(self.index_BPM[j]),int(self.index_BPM[i]),int(self.index_BPM[N])]),2))
                self.alpha_massive2=np.append(self.alpha_massive2,self.beta_from_phase_3(self.coordinats_massive[j],self.coordinats_massive[i],self.coordinats_massive[N],np.array([int(self.index_BPM[j]),int(self.index_BPM[i]),int(self.index_BPM[N])]),12))
        return self.beta_massive1,self.beta_massive0,self.beta_massive2,self.alpha_massive1,self.alpha_massive0,self.alpha_massive2,self.index_BPM,self.freq
