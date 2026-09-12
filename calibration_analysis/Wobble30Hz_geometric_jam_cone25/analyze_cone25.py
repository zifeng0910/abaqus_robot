from pathlib import Path
src=Path(__file__).resolve().parents[1]/'Wobble30Hz_wallon_8p333'/'analyze_wobble_8p333.py'
txt=src.read_text(encoding='utf-8')
txt=txt.replace("JOB='Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival'; OUT=HERE", "JOB='Wobble_F30_G6L45_WallOn_Free_0083_Cone25_WobbleSurvival'; OUT=Path(__file__).resolve().parent")
exec(compile(txt,str(src),'exec'),globals(),globals())
