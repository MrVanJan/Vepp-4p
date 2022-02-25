import numpy as np


class Calculate_Freq():
    '''
    класс Calculate_Freq - возвращает частоту соответствующую максимуму спектральной плотности сигнала
    x-входной массив значений координат
    в качестве оконной функции используется окно ханинга
    интерполяция спектра производится пораболой по 3-м точкам
    для этого метод find_local_maxes - находит фурье и максимум спектра,после этого методом calculate_freq_parabula_interpolation
    производится интерполяция спектра.
    '''
    def __init__(self,x:np.array):
        self.x=x

    @property
    def find_local_maxes(self):
        self.maxes_array = np.empty(0)
        self.WF = np.hanning(self.x.shape[0])
        self.freqArray = np.abs(np.fft.fft(self.x*self.WF))
        for i in np.arange(1,self.x.shape[0]-1).tolist():
            if (self.freqArray[i] > self.freqArray[i - 1]) and (self.freqArray[i] > self.freqArray[i + 1]):
                self.maxes_array= np.append(self.maxes_array,[i, self.freqArray[i]])

        self.maxes_array=np.reshape(self.maxes_array,(int(self.maxes_array.shape[0]/2),2))
        self.maxes_array = self.maxes_array[np.argsort(self.maxes_array[:, 1])]
        self.major_max=int(self.maxes_array[-1][0])
        return self.major_max

    def calculate_freq_parabola_interpolation(self):
        '''
        Решаю систему линейных уравнений типа aXi^2+bXi+c=Yi
        A.(a,b,c)=(Y0,Y1,Y2)
        '''
        self.A_matrix=np.array([[(self.major_max-1)**2,(self.major_max-1),1]
                                        ,[self.major_max**2,self.major_max,1]
                                        ,[(self.major_max+1)**2,(self.major_max+1),1]])
        self.Y_matrix=np.array([self.freqArray[self.major_max-1],self.freqArray[self.major_max],self.freqArray[self.major_max+1]])
        self.cooficiets=np.linalg.solve(self.A_matrix,self.Y_matrix)
        self.freq=-self.cooficiets[1]/(2*self.cooficiets[0]*self.freqArray.shape[0])
        if self.freq<0.5:
            self.freq=1-self.freq

        return self.freq

    @property
    def calculate_freq(self):
        self.find_local_maxes
        return self.calculate_freq_parabola_interpolation()


    def calculate_ampl_and_phase(self,len):
        self.freq=self.calculate_freq
        #self.freq=0.5368830987373664
        self.sumWF=np.sum(np.hanning(len))
        self.C=(2/self.sumWF)*np.sum(np.hanning(len)*self.x[:len]*np.cos([2*np.pi*self.freq*i for i in np.arange(len).tolist()]))
        self.S=(2/self.sumWF)*np.sum(np.hanning(len)*self.x[:len]*np.sin([2*np.pi*self.freq*i for i in np.arange(len).tolist()]))
        self.ampl=np.sqrt(self.C**2+self.S**2)
        self.phase=np.arctan2(-self.S, self.C)
        return self.ampl,self.phase