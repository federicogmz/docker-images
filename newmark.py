#!/usr/bin/env python
#
##############################################################################
#
# MODULE:       Newmark
#
# AUTHOR(S):    Federico Gomez
#
# PURPOSE:      Newmark displacements (Jibson, 2007)
#
# DATE:         Wed Oct  4 10:02:02 2017
#
##############################################################################

#%module
#% description: Newmark displacements
#%end

#%option G_OPT_R_INPUT
#% key: dem
#% key_desc: dem
#% required: yes
#% description: Digital elevation model
#%end
#%option G_OPT_V_INPUT
#% key: mask
#% key_desc: mask
#% required: yes
#% description: Mask
#%end

#%option G_OPT_R_INPUT
#% key: c
#% key_desc: c
#% required: yes
#% description: Effective cohesion ()
#%end
#%option G_OPT_R_INPUT
#% key: gamma
#% key_desc: gamma
#% required: yes
#% description: Material unit weight ()
#%end
#%option G_OPT_R_INPUT
#% key: z
#% key_desc: z
#% required: yes
#% description: Soil thickness (m)
#%end
#%option G_OPT_R_INPUT
#% key: phy
#% key_desc: phy
#% required: yes
#% description: Effective friction angle ()
#%end

#%option G_OPT_R_INPUT
#% key: cf_cohe
#% description: Cohesion correction factor ()
#% type: double
#% required: yes
#% key_desc: cf_cohe
#% answer: 1
#%end
#%option G_OPT_R_INPUT
#% key: cf_gamma
#% description: Material unit weight correction factor ()
#% type: double
#% required: yes
#% key_desc: cf_gamma
#% answer: 1
#%end
#%option G_OPT_R_INPUT
#% key: cf_z
#% description: Soil thickness correction factor()
#% type: double
#% required: yes
#% key_desc: cf_z
#% answer: 1
#%end
#%option G_OPT_R_INPUT
#% key: cf_phy
#% description: Friction correction factor ()
#% type: double
#% required: yes
#% key_desc: cf_phy
#% answer: 1
#%end

#%option G_OPT_R_INPUT
#% key: m
#% description: Proportion of the slab thickness that is saturated (0-1)
#% type: double
#% required: yes
#% key_desc: m
#% answer: 1
#%end
#%option G_OPT_R_INPUT
#% key: pga
#% description: Peak ground acceleration
#% type: double
#% required: yes
#% key_desc: pga
#%end

#%option G_OPT_R_OUTPUT
#% description: Newmark displacements (cm)
#% key: nd
#% gisprompt: new, cell, raster
#% required: yes
#% answer: newmark
#%end

#import modules to be used
import sys
import atexit
import grass.script as grass
from grass.script import parser, run_command

#define function to clean up intermediate maps
def cleanup():
    run_command('g.remove', flags='f', type='raster',
                name=('slope_deg','t','cohesive','frictional','pp',
                      'fs','ac','nd','ndc','MASK'))

def main():
    #Set computationa region from mask 
    run_command("g.region",
                vector=options['mask'])
    run_command("r.mask",
                vector=options['mask'])

    #Calculate slope map from dem
    run_command("r.slope.aspect",
                elevation = options['dem'],
                slope = "slope_deg",
                format = "degrees",
                precision = "FCELL",
                zscale = 1.0,
                min_slope = 0.0)

    #Calculate soil thickness normal to slope
    grass.mapcalc("t=$cf_z*$z*cos(slope_deg)",
                  z = options['z'],
                  cf_z = options['cf_z'])

    #Calculate cohesive factor of fs
    grass.mapcalc("cohesive=$c*$cf_cohe/($cf_gamma*$gamma*t*sin(slope_deg))",
                  c = options['c'],
                  cf_cohe = options['cf_cohe'],
                  gamma = options['gamma'],
                  cf_gamma = options['cf_gamma'])

    #Calculate frictional factor of fs
    grass.mapcalc("frictional=tan($cf_phy*$phy)/tan(slope_deg)",
                  phy = options['phy'],
                  cf_phy = options['cf_phy'])

    #Calculate pore pressure factor of fs
    grass.mapcalc("pp=($m*$gammaw*tan($cf_phy*$phy))/($cf_gamma*$gamma*tan(slope_deg))",
                  m = options['m'],
                  phy = options['phy'],
                  cf_phy = options['cf_phy'],
                  gamma = options['gamma'],
                  cf_gamma = options['cf_gamma'],
                  gammaw = 9.81)

    #Calculate Factor of Safety
    grass.mapcalc("fs=cohesive+frictional-pp")

    #Calculate Critical Aceleration
    grass.mapcalc("ac=(fs-1)*sin(slope_deg)")

    #Calculate Newmark displacements
    grass.mapcalc("nd= exp(10,(0.215+ log((exp((1-(ac/$pga)),2.341)*exp((ac/$pga),-1.438)),10)))",
                  pga=options['pga'])

    #Reclassify newmark displacements
    grass.mapcalc("ndc=if(nd>=0 && nd<1,2,if(nd>=1 && nd<5,3,if(nd>=5 && nd<15,4,5)))")

    #Join to get output map
    grass.mapcalc("newmark=if((fs<1),6,if((ac>$pga),1,ndc))",
                  pga = options['pga'])

    #Join to get output map
    grass.mapcalc("nw_prob=0.335*(1-exp((-0.048*exp(nd,1.565))))")

    #Set colors
    colors_rules = "6 255:0:0\n5 252:127:0\n4 255:255:0\n3 127:255:0\n2 5:157:2\n1 10:60:5"
    grass.write_command('r.colors', map='newmark', rules='-',
                        stdin=colors_rules, quiet=True)

    #Rename output map
    output=options['nd']
    run_command('g.rename',
                raster = ("newmark", output))

    return 0

if __name__ == "__main__":
    options, flags = parser()
    atexit.register(cleanup)
    sys.exit(main())
