import numpy as np


class Method_two_pickup():
    '''
    класс Method_two_pickup
    x_massive1,x_massive2 - массивы координат с двух различных датчиков
    beta1,beta2,alph1 - значения опт.ф-ий на двух различных датчиках
    delt_phase - набек между двумя датчиками

    self.x_massive1[0].shape[0] - число оборотов
    liniarization_trajectory - метод линеаризации траектории,не стал делать как отдельный класс т.к линеаризация нужна только после нахождения не линеаризованной траектории.
    '''

    def __init__(self,x_massive1:np.array,x_massive2:np.array,beta1:float,beta2:float,delt_phase,alph1:float):
        self.x_massive1=x_massive1
        self.x_massive2 = x_massive2
        self.beta1=beta1
        self.beta2 = beta2
        self.delt_phase=delt_phase
        self.alph1=alph1
    def px_massive(self):
        self.coordinats_Array=np.empty(0)
        for j in np.arange(self.x_massive1.shape[0]).tolist():
            self.x=self.x_massive1[j]
            self.p=(self.x_massive2[j]-np.sqrt(self.beta2/self.beta1)*(np.cos(self.delt_phase)+self.alph1*np.sin(self.delt_phase)))/(np.sqrt(self.beta2*self.beta1)*np.sin(self.delt_phase))
            self.coordinats_Array=np.append(self.coordinats_Array,np.array([self.x,self.p]))
        self.vector_massive = np.reshape(self.coordinats_Array, (self.x_massive1.shape[0], 2))
        return np.reshape(self.coordinats_Array,(self.x_massive1.shape[0],2)).T

    def liniarization_trajectory(self,beta1_model,alf1_model):
        self.normalization_matrix = np.array([[1 / np.sqrt(beta1_model), 0],[alf1_model/ np.sqrt(beta1_model),np.sqrt(beta1_model)]])
        self.norm_coordinats_Array = np.dot(self.normalization_matrix, self.vector_massive.T).T
        return np.reshape(self.norm_coordinats_Array, (self.x_massive1.shape[0], 2)).T