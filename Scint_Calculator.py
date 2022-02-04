#WebPage code to calculate Scint bandwidth and timescale
#Made by Juan G. Lebron Medina
# import libraries
from flask import Flask,render_template,request
import astropy.units as u
import astropy.coordinates as c
import pygedm
#Web page code
app= Flask(__name__)
@app.route("/")
def main():
    return render_template("compute.html")
@app.route("/send", methods=["POST"])
def send():
    #Getting values from page
    if request.method == "POST":
        RA = (request.form["RA"])
        DEC = (request.form["DEC"])
        DM = float(request.form["DM"])
        nu = (request.form["nu"])
        V_perp = (request.form["V"])
        #In case user dont input nu or V
        if nu == "":
            nu=1
        else: 
            nu = float(nu)
        if V_perp == "":
            V_perp=100
        else: 
            V_perp= float(V_perp)



        # Code to calculate scint
        ############################
        # coordinates and distance #
        ############################
        # add units 
        ra = c.Angle(RA, unit='hourangle') # right ascension
        dec = c.Angle(DEC, unit='degree') # declination
        # create sky coordinates
        sky_coords = c.SkyCoord(ra, dec, frame='icrs')

        #Pygedm version 3.3.0
        Array= pygedm.ne2001_wrapper.dm_to_dist(sky_coords.galactic.l.value,
                                        sky_coords.galactic.b.value,
                                        DM,nu,
                                        full_output=True)

        #find smtau, distance, scattering time
        smtau = (Array.get('smtau'))
        d = u.Quantity((Array.get('dist')), unit='kpc')
        tauiss = u.Quantity((Array.get('tau_sc')), unit='s')
        sm=smtau #tau uses difrent units 


        #d in (pc) and tauiss in (s)
        # convert distance to kpc
        d_kpc = d.to('kpc')
        d_kpc = round(d_kpc.value,2) #Lazy rounding
        # convert scattering timescale to ms
        tauiss_ms = tauiss.to('ms')
        
        #Max DM distance multiply by 100
        LIST = pygedm.ne2001_wrapper.dist_to_dm(sky_coords.galactic.l.value,
                                                sky_coords.galactic.b.value,
                                                100*d.value,nu,
                                                full_output=True)

        Max_DM = round((LIST.get('dm')),2)

        ###########################
        # scintillation bandwidth #
        ###########################

        c1 = 1.16 # param for uniform, Kolmogorov medium

        # calculate scintillation bandwidth (scintbw, kHz)
        # Ref: Eq. 35 from Cordes & Lazio 1991
        scintbw = u.Quantity(c1/(2*3.1415926535*tauiss_ms.value), unit='kHz')

        ############################
        # scintillation time scale #
        ############################
        # user-defined quantities

        # pulsar perpendicular/transverse velocity
        Units = u.Quantity(3.3, unit='s')
        # calculate scintillation time scale (scinttime, s)
        # Ref: Eq. 46 from Cordes & Lazio 1991 (2.3 change to 3.3)
        scinttime = u.Quantity(Units*nu**1.2*sm**(-0.6)*(100/V_perp), unit='s')
        scintbw = round(scintbw.value,1)
        scinttime = round(scinttime.value,1)



        return render_template("compute.html",ra=ra, dec=dec,DM=DM,nu=nu, V_perp=V_perp, d_kpc=d_kpc, Max_DM=Max_DM,  scintbw=scintbw, scinttime= scinttime)




if __name__ == "__main__":
    app.run()
    
