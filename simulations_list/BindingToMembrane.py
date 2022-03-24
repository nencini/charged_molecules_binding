import numpy as np
import sys
import gc
import math
import os
import shutil

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import


from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
plt.rcParams["figure.figsize"] = (10,10)








class AnalyzeBindingDefinition(object):
    def __init__(self,road,name,index,particles,size,limits,beginTime):
        self.road = road
        self.name = name
        self.index = index
        self.particles = particles
        self.size=size
        self.limits=limits
        self.beginTime=beginTime
        should_continue=1


        if not os.path.isdir('./'+str(name)):
            os.mkdir(name)

        if self.index=="yes":
            self.create_index()

        if self.index=="SMS":
            should_continue=self.create_index_SMS()

        if should_continue==1:
            for cutOff in np.arange(limits[0],limits[1],limits[2]):
                gc.collect()
                self.calculate_binding(cutOff)


    def create_index(self):
        os.system('echo "del 0-1 \n del 2-3 \n splitres 1 \n del 1 \n q \n \n "|gmx make_ndx -f ' + self.road + '/' + self.name + '/' + self.name + '.tpr -o ' + self.name + '/' + self.name + '.ndx')


    def create_index_SMS(self):
        os.system('echo "q \n \n "|gmx make_ndx -f ' + self.road + '/' + self.name + '/' + self.name + '.tpr -o ' + self.name + '/aux.ndx')
        should_continue=input("Is the following satisfied? \n 2 POPC \n 3,4,5 non-standard residues \n 8 Protein \n \n Confirm by YES ")

        if should_continue=="YES":
            os.system('echo "del 0-1 \n del 4-5\n del 5-14 \n 1|2|3|4 \n splitch 5 \n del 1-5 \n q \n \n "|gmx make_ndx -f ' + self.road + '/' + self.name + '/' + self.name + '.tpr -o ' + self.name + '/' + self.name + '.ndx')

            satisfied=input("Did the spliting succeeded? YES/NO")
            if satisfied=="YES":
                return 1
            else:
                return 0
        else:
            print("Conditions are not met")
            return 0



    def calculate_binding(self,cutOff):
        cutOut=int(round(1000*cutOff))
        forRM=str(self.name)+"_"+str(cutOut)+".dat"
        if os.path.exists(forRM):
            os.remove(forRM)
        os.system('echo "0 \n 1 \n 2 \n 3 \n 4 \n 5 \n 6 \n 7 \n 8 \n 9 \n 10 \n 11 \n 12 \n 13 \n 14 \n 15 \n 16 \n 17 \n 18 \n 19 \n 20 \n 21\n 22\n 23 \n 24\n 25 \n 26 \n 27 \n 28 \n 29 \n 30 \n 31\n 32\n 33 \n 34\n 35 \n 36 \n 37 \n 38 \n 39 \n 40 \n 41\n 42\n 43 \n 44\n 45 \n 46 \n 47 \n 48 \n 49 \n 50 \n 51\n 52 \n 53 \n 54\n 55 \n 56 \n 57 \n 58 \n 59 \n 60 \n 61 \n 62 \n 63 \n 64\n 65 \n 66 \n 67 \n 68 \n 69 \n 70 \n 71\n 72 \n 73 \n 74\n 75 \n 76 \n 77 \n 78 \n 79 \n 80 \n 81\n 82 \n 83 \n 84\n 85 \n 86 \n 87 \n 88 \n 89 \n 90 \n 91\n 92 \n 93 \n 94\n 95 \n 96 \n 97 \n 98 \n 99 \n 100 \n 101\n 102 \n 103 \n 104\n 105 \n 106 \n 107 \n 108 \n 109 \n 110 \n 111\n 112 \n 113 \n 114\n 115 \n 116 \n 117 \n 118 \n 119 \n 120 \n 121\n 122 \n 123 \n 124\n 125 \n 126 \n 127 \n 128 \n 129 \n 130 \n 131\n 132 \n 133 \n 134\n 135 \n 136 \n 137 \n 138 \n 139 \n 140 \n 141\n 142 \n 143 \n 144\n 145 \n 146 \n 147 \n 148 \n 149 \n 150  \n 151 \n 152 \n 153 \n 154 \n 155 \n 156 \n 157 \n 158 \n 159 \n 160 \n 161 \n 162 \n 163 \n 164 \n 165 \n 166 \n 167 \n 168 \n 169 \n 170 \n 171\n 172\n 173 \n 174\n 175 \n 176 \n 177 \n 178 \n 179 \n 180 \n 181\n 182\n 183 \n 184\n 185 \n 186 \n 187 \n 188 \n 189 \n 190 \n 191\n 192\n 193 \n 194\n 195 \n 196 \n 197 \n 198 \n 199 \n 200 \n 201\n 202 \n 203 \n 204\n 205 \n 206 \n 207 \n 208 \n 209 \n 210 \n 211 \n 212 \n 213 \n 214\n 215 \n 216 \n 217 \n 218 \n 219 \n 220 \n 221\n 222 \n 223 \n 224\n 225 \n 226 \n 227 \n 228 \n 229 \n 230 \n 231\n 232 \n 233 \n 234\n 235 \n 236 \n 237 \n 238 \n 239 \n 240 \n 241\n 242 \n 243 \n 244\n 245 \n 246 \n 247 \n 248 \n 249 \n 250 \n 251\n 252 \n 253 \n 254\n 255 \n 256 \n 257 \n 258 \n 259 \n 260 \n 261\n 262 \n 263 \n 264\n 265 \n 266 \n 267 \n 268 \n 269 \n 270 \n 271\n 272 \n 273 \n 274\n 275 \n 276 \n 277 \n 278 \n 279 \n 280 \n 281\n 282 \n 283 \n 284\n 285 \n 286 \n 287 \n 288 \n 289 \n 290 \n 291\n 292 \n 293 \n 294\n 295 \n 296 \n 297 \n 298 \n 299 \n 300 \n 301 \n 302 \n 303 \n 304 \n 305 \n 306 \n 307 \n 308 \n 309 \n 310 \n 311 \n 312 \n 313 \n 314 \n 315 \n 316 \n 317 \n 318 \n 319 \n 320 \n 321\n 322\n 323 \n 324\n 325 \n 326 \n 327 \n 328 \n 329 \n 330 \n 331\n 332\n 333 \n 334\n 335 \n 336 \n 337 \n 338 \n 339 \n 340 \n 341\n 342\n 343 \n 344\n 345 \n 346 \n 347 \n 348 \n 349 \n 350 \n 351\n 352 \n 353 \n 354\n 355 \n 356 \n 357 \n 358 \n 359 \n 360 \n 361 \n 362 \n 363 \n 364\n 365 \n 366 \n 367 \n 368 \n 369 \n 370 \n 371\n 372 \n 373 \n 374\n 375 \n 376 \n 377 \n 378 \n 379 \n 380 \n 381\n 382 \n 383 \n 384\n 385 \n 386 \n 387 \n 388 \n 389 \n 390 \n 391\n 392 \n 393 \n 394\n 395 \n 396 \n 397 \n 398 \n 399 \n 400 \n 401\n 402 \n 403 \n 404\n 405 \n 406 \n 407 \n 408 \n 409 \n 410 \n 411\n 412 \n 413 \n 414\n 415 \n 416 \n 417 \n 418 \n 419 \n 420 \n 421\n 422 \n 423 \n 424\n 425 \n 426 \n 427 \n 428 \n 429 \n 430 \n 431\n 432 \n 433 \n 434\n 435 \n 436 \n 437 \n 438 \n 439 \n 440 \n 441\n 442 \n 443 \n 444\n 445 \n 446 \n 447 \n 448 \n 449 \n 450"|gmx mindist -f ' + self.road + '/' + self.name + '/' +self.name + '.xtc -s ' + self.road + '/' + self.name + '/' + self.name + '.tpr -ng ' + str(self.particles) + ' -on ' + self.name + '/' + self.name + '_' + str(cutOut) + ' -n ' + self.name + '/' + self.name +'.ndx -d ' + str(cutOff) + ' -group -b ' + str(self.beginTime))

        #if os.path.exists(str(self.name)+"_"+str(cutOut)+".xvg"):
        #    os.rename(str(self.name)+"_"+str(cutOut)+".xvg",str(self.name)+"/"+str(self.name)+"_"+str(cutOut)+".xvg")

        if os.path.exists(str(self.name)+"/"+str(self.name)+"_"+str(cutOut)+".dat"):
            os.remove(str(self.name)+"/"+str(self.name)+"_"+str(cutOut)+".dat")
        for x in range (1,self.size+1):
            self.try_different_number_of_atoms(x,cutOut)

        if os.path.exists("mindist.xvg"):
            os.remove("mindist.xvg")



    def try_different_number_of_atoms(self,number,cutOut):
        pure=np.loadtxt(str(self.name)+"/"+str(self.name)+"_"+str(cutOut)+".xvg", comments=["@","#"])
        result = (pure[:,1:] > number).astype(int)
        bound=np.count_nonzero(result, axis=1)
        np.savetxt(str(self.name)+"_"+str(cutOut)+"_"+str(number)+".txt",bound)
        forAwk=str(self.name)+"_"+str(cutOut)+"_"+str(number)+".txt"

        average=0
        with open(forAwk) as fpin:
            for nr, line in enumerate(fpin):
                numerical=float(line)
                average+=numerical
                last=nr


        with open(str(self.name)+"/"+str(self.name)+"_"+str(cutOut)+".dat", 'a') as f:
            f.write("%s %.4f \n" % (number, average/last))

        os.remove(forAwk)


class Time_evolution(object):
    #there are 3 different modes one can choose from: evolve, cut, inflex
    #evolve - calculates binding in time for given cut-offs and for all posible # of atoms in molecule
    #cut - for the given range of cut-offs, it creates files for each posible # of atoms and writes down the bindings in different cut-offs
    def __init__(self,name,size,what,limits):
        self.size=size
        self.name = name
        self.what=what
        self.limits=limits

        if what=="evolve":
            if not os.path.isdir('./'+str(name)+'/inTime'):
                os.mkdir(str(name)+'/inTime')


        if what=="cut":
            if not os.path.isdir('./'+str(name)+'/cut'):
                os.mkdir(str(name)+'/cut')
            else:
                shutil.rmtree(str(name)+'/cut')
                os.mkdir(str(name)+'/cut')

        if what=="inflex":
            if not os.path.isdir('./'+str(name)+'/derivative'):
                os.mkdir(str(name)+'/derivative')
            if os.path.exists(str(self.name)+'/derivative/minima.dat'):
                os.remove(str(self.name)+'/derivative/minima.dat')

        if 1==1:
            for cutOff in np.arange(limits[0],limits[1],limits[2]):
                self.calculate_evol(cutOff)


        if 1==2:
            for cutOff in np.arange(0.2,0.275,0.025):
                self.calculate_evol(cutOff)


            self.calculate_evol(0.27)
            self.calculate_evol(0.275)
            self.calculate_evol(0.28)
            self.calculate_evol(0.29)
            self.calculate_evol(0.30)
            self.calculate_evol(0.31)
            self.calculate_evol(0.32)
            self.calculate_evol(0.325)
            self.calculate_evol(0.33)
            self.calculate_evol(0.34)
            self.calculate_evol(0.35)
            self.calculate_evol(0.36)
            self.calculate_evol(0.37)
            self.calculate_evol(0.375)
            self.calculate_evol(0.38)
            self.calculate_evol(0.39)

            for cutOff in np.arange(0.4,0.7,0.025):
                self.calculate_evol(cutOff)



    def calculate_evol(self,cutOff):
        if self.what=="evolve":
            for number in range (1,self.size+1):
                cutOut=int(round(1000*cutOff))
                pure=np.loadtxt(str(self.name)+"/"+str(self.name)+"_"+str(cutOut)+".xvg", comments=["@","#"])
                result = (pure[:,1:] > number).astype(int)
                bound=np.count_nonzero(result, axis=1)
                np.savetxt(str(self.name)+'/inTime/'+str(self.name)+"_"+str(cutOut)+"_"+str(number)+".txt",bound)

        if self.what=="cut":
            for number in range (0,self.size,1):
                cutOut=int(round(1000*cutOff))
                pure=np.loadtxt(str(self.name)+"/"+str(self.name)+"_"+str(cutOut)+".dat", comments=["@","#"])
                result = pure[number][1]
                with open(str(self.name)+'/cut/'+str(self.name)+"_"+str(number)+".dat", 'a') as f:
                    f.write("%.3f %.4f \n" % (cutOff, result))

        if self.what=="inflex":
            cutOut=int(round(1000*cutOff))
            pure=np.loadtxt(str(self.name)+"/"+str(self.name)+"_"+str(cutOut)+".dat", comments=["@","#"])
            derivative=pure.copy()
            secondDerivative=pure.copy()
            i=0
            for number in range (0,self.size-1,1):
                derivative[number][1]=pure[number+1][1]-pure[number][1]
                if derivative[number][1]<i:
                    i=derivative[number][1]
                    index=number

            print(pure[index][1])

            with open(str(self.name)+'/derivative/minima.dat', 'a') as f:
                    f.write("%.3f %.4f \n" % (cutOff, pure[index][1]))
            np.savetxt(str(self.name)+'/derivative/'+str(self.name)+"_"+str(cutOut)+"_deriv.dat",derivative)

            for number in range (0,self.size-1,1):
                secondDerivative[number][1]=derivative[number+1][1]-derivative[number][1]
            np.savetxt(str(self.name)+'/derivative/'+str(self.name)+"_"+str(cutOut)+".dat",secondDerivative)


class PlotOPvsbinding(object):
    def __init__(self,name,mol):
        self.name = name
        binding=np.loadtxt(str(self.name)+'/'+str(self.name)+"_bining.dat",usecols=[13])
        op= np.loadtxt('../OP/final_dat_data/'+str(self.name)+".dat",usecols=[1,5,7])
        with open('../relationship/'+str(mol)+".dat", 'a') as f:
            f.write("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f \n" % (binding[0], np.average(binding),np.std(binding),op[np.where(op[:,0] == 1)[0][0],1],op[np.where(op[:,0] == 1)[0][0],2],op[np.where(op[:,0] == 2)[0][0],1],op[np.where(op[:,0] == 2)[0][0],2],op[np.where(op[:,0] == 3)[0][0],1],op[np.where(op[:,0] == 3)[0][0],2],op[np.where(op[:,0] == 4)[0][0],1],op[np.where(op[:,0] == 4)[0][0],2]))

