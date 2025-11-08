
import numpy as np
from scipy.constants import e
from scipy.optimize import curve_fit

class HRMR_PPMS():
    def __init__(self, path: str, channel=[1,2]):
        assert path[-4:] == ".dat", "Input files must be .dat extension"
        self.path = path
        skiprows = 0
        with open(path, "r") as f:
            for i in range(70):
                if f.readline() == "[Data]\n":
                    skiprows = i + 2
                    break
            else:
                raise ValueError("The measured data was not found in the input file.")
        assert channel, "Please specify the channel"
        if type(channel) == int:
            channel = [channel]
        col = [4] + list(map(lambda x:x+18, channel)) # 4: magnetic field (Oe)
        data = np.loadtxt(path , skiprows=skiprows, delimiter=',', usecols=col)
        self.data = data.T
        self.data[0] /= 1e4 # Oe -> T
        self.unit='ohm'

    def remove_even(self, poly=1): # 偶関数の除去。デフォルトは2次関数のみ。poly*2次の関数まで除去。
        coeff_hr = np.polyfit(self.data[0], self.data[2], 2*poly)
        new_data = self.data[2]
        for i in range(poly+1):
            new_data = new_data - self.data[0] ** ((poly - i)*2) * coeff_hr[i*2]
        self.data[2] = new_data
    
    def odd(self):
        self.data[2] = (self.data[2] - np.roll(self.data[2], shift=len(self.data[2])//2)) / 2.

    def to_resistivity(self, thickness_nm, L_per_W = 4):
        self.data[1] = self.data[1] / L_per_W * thickness_nm * 1e-9 * 1e6 # u*Ohm*m
        self.data[2] = self.data[2] * thickness_nm * 1e-9 * 1e6 # u*Ohm*m
        self.unit = "resistivity"

    def mr(self):
        return self.data[1] / self.data[1].min() * 100 - 100

    def _sigma(self): # キャリア密度の計算のためにrho_xx, rho_yxに変換
        assert self.unit=="resistivity", "Please convert the units to resistivity before calculating sigma."
        rho_xx = self.data[1] * 1e-6 # Ohm * m
        rho_yx = self.data[2] * 1e-6 # Ohm * m
        sigma_xx = rho_xx / (rho_xx ** 2 + rho_yx ** 2)
        sigma_yx = - rho_yx / (rho_xx ** 2 + rho_yx ** 2)
        return self.data[0], sigma_xx, sigma_yx
    
    @staticmethod
    def multi_carrier_simu(h, *param): # h1: 磁場(rho_xx), h2: 磁場(rho_yx) paramはキャリア濃度と移動度を交互に入力 (/cm³, cm²/Vs)。pキャリアは両方符号をマイナスにする。
        assert len(param) % 2 == 0, "*param is input carrier concentration and mobility alternately."
        n_carrier = len(param) // 2
        sigma_xx_sigma_yx = np.zeros_like(h)
        for i in range(n_carrier):
            n = param[i * 2 + 0]
            u = param[i * 2 + 1]
            sigma_xx_sigma_yx[:len(h)//2] += n*e*u / (1 + u**2 * h[:len(h)//2]**2)
            sigma_xx_sigma_yx[len(h)//2:] += n*e*u**2*h[len(h)//2:] / (1 + u**2 * h[len(h)//2:]**2)
        return sigma_xx_sigma_yx
    
    @staticmethod
    def multi_carrier_fit(h1, sigma_xx, sigma_yx, p0, **kwargs):
        x = np.concatenate((h1, h1))
        y = np.concatenate((sigma_xx, sigma_yx))
        popt, pcov = curve_fit(HRMR_PPMS.multi_carrier_simu, x, y, p0=p0, **kwargs)
        return popt, pcov
