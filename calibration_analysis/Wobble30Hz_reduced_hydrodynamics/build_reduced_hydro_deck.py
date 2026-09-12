from pathlib import Path
import re
ROOT=Path(r'J:\\abaqusfangzhen'); OUT=Path(__file__).resolve().parent
src=ROOT/'Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival.inp'; job='Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083'
txt=src.read_text(encoding='utf-8',errors='ignore')
# Remove the complete Eulerian fluid part.
txt=re.sub(r'\*Part, name=FLUID_EULERIAN.*?\*End Part\s*', '', txt, flags=re.I|re.S, count=1)
# Remove the fluid instance and all fluid-only assembly sets/initialization.
txt=re.sub(r'\*Instance, name=Fluid_EULERIAN-1.*?\*End Instance\s*', '', txt, flags=re.I|re.S, count=1)
txt=re.sub(r'\*Elset, elset=FLUID_CEL_ALL, instance=Fluid_EULERIAN-1, generate\s*.*?(?=\*Elset|\*Surface|\*End Assembly)', '', txt, flags=re.I|re.S, count=1)
txt=re.sub(r'\*Elset, elset=FLUID_CEL_INIT_ALL, instance=Fluid_EULERIAN-1\s*.*?(?=\*Elset|\*Surface|\*End Assembly)', '', txt, flags=re.I|re.S, count=1)
txt=re.sub(r'\*Surface, type=EULERIAN MATERIAL, name=FLUID_CEL_SURF\s*.*?(?=\*End Assembly)', '', txt, flags=re.I|re.S, count=1)
txt=re.sub(r'\*Initial Conditions, type=VOLUME FRACTION\s*.*?(?=\*|\*MATERIALS)', '', txt, flags=re.I|re.S, count=1)
txt=re.sub(r'\*Material, name=MAT_FLUID_CEL\s*.*?(?=\*Material, name=MAT_PIPE_RIGID)', '', txt, flags=re.I|re.S, count=1)
txt=re.sub(r'^\*\*.*FLUID_EULERIAN.*\n', '', txt, flags=re.I|re.M)
txt=re.sub(r'\*Contact Inclusions\s*([^*]*?)\*Contact Property Assignment\s*([^*]*?)\*', lambda m: '*Contact Inclusions\nROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\n*Contact Property Assignment\nROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF, PROP_CEL_HARD\n*', txt, flags=re.I|re.S, count=1)
txt=txt.replace('** CEL fill uses 1564 verified Eulerian elements; unassigned volume is VOID.','** Reduced-Hydro: CEL fluid part, material, initialization and CEL contact were removed.')
txt=txt.replace('** SOCKET MAGNETIC LOAD INTERFACE -- RP_ROBOT translations and UR1-UR3 free.','** REDUCED-HYDRO LOAD INTERFACE: magnetic Socket loads plus dissipative body-following hydro loads.')
txt=txt.replace('** Six VUAMP amplitudes remain updated every Explicit increment by vuforc_socket_bridge.f.','** SOCKET_* remain production Magpylib loads; HYDRO_* are local dissipative loads from RP V/VR sensors.')
# Add hydro amplitudes and CLOADs immediately after magnetic MZ load.
needle='RP_ROBOT, 6, 1.\n** '
insert='''RP_ROBOT, 6, 1.\n*Amplitude, name=HYDRO_FX, definition=USER\n*Amplitude, name=HYDRO_FY, definition=USER\n*Amplitude, name=HYDRO_FZ, definition=USER\n*Amplitude, name=HYDRO_MX, definition=USER\n*Amplitude, name=HYDRO_MY, definition=USER\n*Amplitude, name=HYDRO_MZ, definition=USER\n*Cload, amplitude=HYDRO_FX\nRP_ROBOT, 1, 1.\n*Cload, amplitude=HYDRO_FY\nRP_ROBOT, 2, 1.\n*Cload, amplitude=HYDRO_FZ\nRP_ROBOT, 3, 1.\n*Cload, amplitude=HYDRO_MX\nRP_ROBOT, 4, 1.\n*Cload, amplitude=HYDRO_MY\nRP_ROBOT, 5, 1.\n*Cload, amplitude=HYDRO_MZ\nRP_ROBOT, 6, 1.\n** '''
txt=txt.replace(needle,insert,1)
# Add velocity/rotation-rate sensors after the existing UR3 sensor.
needle='''*Output, history, frequency=1, sensor, name=RP_UR3\n*Node Output, nset=RP_ROBOT\nUR3\n'''
extra=needle+'''*Output, history, frequency=1, sensor, name=RP_V1\n*Node Output, nset=RP_ROBOT\nV1\n*Output, history, frequency=1, sensor, name=RP_V2\n*Node Output, nset=RP_ROBOT\nV2\n*Output, history, frequency=1, sensor, name=RP_V3\n*Node Output, nset=RP_ROBOT\nV3\n*Output, history, frequency=1, sensor, name=RP_VR1\n*Node Output, nset=RP_ROBOT\nVR1\n*Output, history, frequency=1, sensor, name=RP_VR2\n*Node Output, nset=RP_ROBOT\nVR2\n*Output, history, frequency=1, sensor, name=RP_VR3\n*Node Output, nset=RP_ROBOT\nVR3\n'''
txt=txt.replace(needle,extra,1)
txt=txt.replace('** FREE-TRANSLATION PARAMETERS:', '** REDUCED-HYDRO BASELINE; CEL REMOVED; production external Magpylib retained.\n** FREE-TRANSLATION PARAMETERS:',1)
# With the CEL deformable elements removed, Explicit needs an explicit user
# time increment; retain the requested per-increment Socket update cadence.
txt=txt.replace('*Dynamic, Explicit\n, 0.008333','*Dynamic, Explicit, DIRECT\n1.0e-7, 0.008333',1)
(ROOT/f'{job}.inp').write_text(txt,encoding='utf-8'); print('WROTE',ROOT/f'{job}.inp')
