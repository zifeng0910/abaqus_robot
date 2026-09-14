"""Python/matplotlib rendering for individual GIFs, mosaics, gallery, and contact sheet.

Figure contract: expose rotation/contact/separation independently from axial translation.
Evidence: synchronized global tube and local cross-section views for every case.
Archetype: image plate + quantitative overlays. No score and no automatic winner.
"""

from pathlib import Path
import html
import json
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

HERE=Path(__file__).resolve().parent; SCREEN=HERE.parent; REPO=SCREEN.parents[1]; ROOT=REPO.parent
GIFS=SCREEN/"gifs"; MOSAICS=SCREEN/"family_mosaics"; GALLERY=SCREEN/"gallery"
CURVE=pd.read_csv(ROOT/"CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv")
XYZ=CURVE[["x_mm","y_mm","z_mm"]].to_numpy(float); RP0=np.array([-7.468174204284,-3.676918015967,-9.550745259298])
LOCAL_CURVE=XYZ[np.linalg.norm(XYZ-RP0,axis=1)<3.5]
FRAMES=61; PIPE_RADIUS=.94
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":7,
                     "axes.linewidth":.7,"svg.fonttype":"none","pdf.fonttype":42})


def interp_pose(df):
    q=np.linspace(0,.2,FRAMES); x=df.cycle_fraction.to_numpy(float); out={}
    assert np.all(np.diff(x)>0), "cycle fraction must increase"
    for c in df.columns:
        if c in ("contact_state",): out[c]=df[c].to_numpy()[np.clip(np.searchsorted(x,q),0,len(df)-1)]
        elif c not in ("time_s",): out[c]=np.interp(q,x,df[c].to_numpy(float))
    out["cycle_fraction"]=q; out["time_s"]=q/float(1/(df.time_s.iloc[-1]/.2))
    return out


def case_data(folder):
    item=json.loads((folder/"case_identity.json").read_text()); item.update(json.loads((SCREEN/"metrics"/(item["case_id"]+"_metrics.json")).read_text()))
    pose=pd.read_csv(folder/"pose_light.csv"); p=interp_pose(pose); p["time_s"]=p["cycle_fraction"]/item["frequency_Hz"]
    prox=pd.read_csv(folder/"proximity_light.csv"); p["sector_deg"]=np.interp(p["time_s"],prox.time_s,np.unwrap(np.radians(prox.nearest_wall_sector_deg)))*180/math.pi
    return item,p


def render_case(folder):
    item,p=case_data(folder); out=GIFS/(item["case_id"]+"_motion.gif")
    fig=plt.figure(figsize=(7.2,3.6),dpi=int("133"),facecolor="white")
    ax=fig.add_subplot(1,2,1,projection="3d"); bx=fig.add_subplot(1,2,2)
    mins=LOCAL_CURVE.min(axis=0)-.5; maxs=LOCAL_CURVE.max(axis=0)+.5
    def update(i):
        ax.cla(); bx.cla()
        ax.plot(LOCAL_CURVE[:,0],LOCAL_CURVE[:,1],LOCAL_CURVE[:,2],color="#8d99a6",lw=7,alpha=.18)
        ax.plot(LOCAL_CURVE[:,0],LOCAL_CURVE[:,1],LOCAL_CURVE[:,2],color="#59636e",lw=.8,alpha=.8)
        head=np.array([p[k][i] for k in ("head_x_mm","head_y_mm","head_z_mm")]); tail=np.array([p[k][i] for k in ("tail_x_mm","tail_y_mm","tail_z_mm")]); rp=np.array([p[k][i] for k in ("rp_x_mm","rp_y_mm","rp_z_mm")])
        b=np.array([p[k][i] for k in ("B_x_T","B_y_T","B_z_T")]); b=b/max(np.linalg.norm(b),1e-30)
        ax.plot([head[0],tail[0]],[head[1],tail[1]],[head[2],tail[2]],color="#0072B2",lw=9,solid_capstyle="round")
        ax.scatter(*head,s=22,color="#D55E00",depthshade=False); ax.quiver(*rp,*b,length=.75,color="#009E73",linewidth=1.5)
        ax.set(xlim=(mins[0],maxs[0]),ylim=(mins[1],maxs[1]),zlim=(mins[2],maxs[2]))
        ax.set_box_aspect(maxs-mins); ax.view_init(elev=23,azim=-58); ax.set_axis_off(); ax.set_title("Global tube view",fontsize=8)
        center=np.array([p[k][i] for k in ("center_x_mm","center_y_mm","center_z_mm")]); e1=np.array([p[k][i] for k in ("e1_x","e1_y","e1_z")]); e2=np.array([p[k][i] for k in ("e2_x","e2_y","e2_z")])
        hp=np.array([np.dot(head-center,e1),np.dot(head-center,e2)]); tp=np.array([np.dot(tail-center,e1),np.dot(tail-center,e2)])
        circle=plt.Circle((0,0),PIPE_RADIUS,facecolor="#d9e2e8",edgecolor="#59636e",alpha=.28,lw=1.2); bx.add_patch(circle)
        bx.plot([hp[0],tp[0]],[hp[1],tp[1]],color="#0072B2",lw=8,solid_capstyle="round"); bx.scatter(*hp,s=30,color="#D55E00",zorder=3,label="HEAD"); bx.scatter(*tp,s=22,color="#0072B2",zorder=3,label="TAIL")
        ang=math.radians(p["sector_deg"][i]); bx.scatter(PIPE_RADIUS*math.cos(ang),PIPE_RADIUS*math.sin(ang),marker="x",s=35,color="#CC79A7",zorder=4)
        bp=np.array([np.dot(b,e1),np.dot(b,e2)]); bx.arrow(0,0,.55*bp[0],.55*bp[1],width=.012,head_width=.08,color="#009E73",length_includes_head=True)
        bx.axhline(0,color="#c5cbd0",lw=.5); bx.axvline(0,color="#c5cbd0",lw=.5); bx.set_aspect("equal"); bx.set(xlim=(-1.08,1.08),ylim=(-1.08,1.08),xlabel="local e1 (mm)",ylabel="local e2 (mm)",title="Local pipe cross-section")
        bx.grid(False)
        fig.suptitle("%s   t=%.3f ms   cycle=%.3f"%(item["case_id"],p["time_s"][i]*1e3,p["cycle_fraction"][i]),fontsize=9,y=.985)
        fig.text(.5,.015,"L=%.2f mm  B=%.0f mT  G=%.0f mT  cone=%.0f deg  f=%.0f Hz  mu=%.2f  zeta=%.2f | tilt=%.1f deg  phase=%+.1f deg  ds=%+.2f mm  %s"%
                 (item["length_mm"],item["B0_mT"],item["gradient_mT"],item["cone_deg"],item["frequency_Hz"],item["mu"],item["zeta"],p["directed_tilt_deg"][i],p["robot_phase_deg"][i],p["delta_s_mm"][i],p["contact_state"][i]),ha="center",va="bottom",fontsize=6.2)
        return []
    animation=FuncAnimation(fig,update,frames=FRAMES,interval=70,blit=False); animation.save(out,writer=PillowWriter(fps=14)); plt.close(fig); print("GIF",out.name,flush=True)


FAMILIES=[
 ("01_geometry_family.gif",[("GEO_200","GEO_200"),("GEO_210","GEO_210"),("GEO_220","GEO_220"),("GEO_230","GEO_230"),("GEO_240","GEO_240")]),
 ("02_gradient_family.gif",[("GRAD_0","GRAD_0"),("GRAD_2","GRAD_2"),("GRAD_4","GRAD_4"),("GRAD_6","GEO_230"),("GRAD_8","GRAD_8")]),
 ("03_cone_family.gif",[("CONE_20","CONE_20"),("CONE_30","GEO_230"),("CONE_40","CONE_40"),("CONE_50","CONE_50")]),
 ("04_B_family.gif",[("B_6","B_6"),("B_10","GEO_230"),("B_14","B_14")]),
 ("05_frequency_family.gif",[("FREQ_20","FREQ_20"),("FREQ_30","GEO_230"),("FREQ_40","FREQ_40")]),
 ("06_friction_family.gif",[("MU_000","MU_000"),("MU_030","GEO_230"),("MU_080","MU_080"),("MU_150","MU_150")]),
 ("07_collision_damping_family.gif",[("ZETA_015","ZETA_015"),("ZETA_030","ZETA_030"),("ZETA_050","GEO_230"),("ZETA_070","ZETA_070")]),
 ("08_hydro_family.gif",[("CPAR_X025","CPAR_X025"),("CPAR_X1","GEO_230"),("CPAR_X4","CPAR_X4"),("KWOB_X025","KWOB_X025"),("KWOB_X1","GEO_230"),("KWOB_X4","KWOB_X4")]),
 ("09_combo_family.gif",[("COMBO_A","COMBO_A"),("COMBO_B","COMBO_B"),("COMBO_C","COMBO_C"),("COMBO_D","COMBO_D"),("COMBO_E","COMBO_E"),("COMBO_F","COMBO_F")])]


def mosaics():
    font=ImageFont.load_default()
    for filename,entries in FAMILIES:
        sources=[Image.open(GIFS/(case+"_motion.gif")) for _,case in entries]; cols=3; rows=math.ceil(len(entries)/cols); frames=[]
        for i in range(FRAMES):
            canvas=Image.new("RGB",(cols*400,rows*220),"white"); draw=ImageDraw.Draw(canvas)
            for j,((label,_),src) in enumerate(zip(entries,sources)):
                src.seek(min(i,src.n_frames-1)); im=src.convert("RGB"); im.thumbnail((400,200)); x=(j%cols)*400; y=(j//cols)*220+20
                canvas.paste(im,(x,y)); draw.text((x+6,(j//cols)*220+4),label,fill="black",font=font)
            frames.append(canvas)
        frames[0].save(MOSAICS/filename,save_all=True,append_images=frames[1:],duration=72,loop=0,optimize=True)
        for src in sources: src.close()
        print("MOSAIC",filename,flush=True)


def contact_sheet(metrics):
    ids=metrics.case_id.tolist(); fig,axes=plt.subplots(math.ceil(len(ids)/4),4,figsize=(16,18),dpi=300); axes=np.asarray(axes).ravel()
    for ax,case_id in zip(axes,ids):
        gif=Image.open(GIFS/(case_id+"_motion.gif")); ims=[]
        for frame in (0,(gif.n_frames-1)//2,gif.n_frames-1): gif.seek(frame); ims.append(np.asarray(gif.convert("RGB").resize((320,160))))
        ax.imshow(np.concatenate(ims,axis=1)); ax.set_title(case_id,fontsize=7,loc="left"); ax.axis("off"); gif.close()
    for ax in axes[len(ids):]:ax.axis("off")
    fig.suptitle("Coarse motion-mode screen: early / mid / late",fontsize=11); fig.tight_layout()
    fig.savefig(SCREEN/"motion_screening_contact_sheet.png",dpi=300,bbox_inches="tight")
    fig.savefig(SCREEN/"motion_screening_contact_sheet.svg",bbox_inches="tight")
    fig.savefig(SCREEN/"motion_screening_contact_sheet.pdf",bbox_inches="tight")
    plt.close(fig)


def gallery(metrics):
    columns=["length_mm","B0_mT","gradient_mT","cone_deg","frequency_Hz","mu","zeta","robot_phase_advance_deg","delta_s_mm","contact_count","longest_contact_us","rotation_translation_index"]
    blocks=[]
    for family,group in metrics.groupby("family",sort=False):
        cards=[]
        for row in group.itertuples():
            table="".join("<tr><th>%s</th><td>%s</td></tr>"%(html.escape(c),html.escape("%.4g"%getattr(row,c) if isinstance(getattr(row,c),(float,np.floating)) else str(getattr(row,c)))) for c in columns)
            cards.append("<article><h3>%s</h3><img loading='lazy' src='../gifs/%s_motion.gif'><table>%s</table></article>"%(row.case_id,row.case_id,table))
        blocks.append("<section><h2>%s</h2><div class='grid'>%s</div></section>"%(family,"".join(cards)))
    page="""<!doctype html><meta charset='utf-8'><title>Coarse motion screen</title><style>body{font:14px Arial,sans-serif;margin:24px;color:#20252a}h1{font-size:22px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:16px}article{border:1px solid #ccd2d8;padding:10px;border-radius:6px}img{width:100%;height:auto}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border-top:1px solid #e4e7ea;padding:3px;text-align:left}th{font-weight:500;color:#59636e}</style><h1>Human-in-the-loop coarse motion-mode screen</h1><p>No score or automatic winner. Review rotation, impact, separation, sliding, tumble, and jam directly.</p>"""+"".join(blocks)
    (GALLERY/"motion_screening_gallery.html").write_text(page,encoding="utf-8")


def review_template(metrics):
    columns=["case_id","visual_rotation_quality","wall_impact_quality","separation_quality","too_much_translation","too_much_sliding","tumble","jam","overall_visual_rank","user_notes"]
    pd.DataFrame([{c:(row.case_id if c=="case_id" else "") for c in columns} for row in metrics.itertuples()]).to_csv(SCREEN/"human_review_template.csv",index=False)


def main():
    for path in (GIFS,MOSAICS,GALLERY):path.mkdir(parents=True,exist_ok=True)
    metrics=pd.read_csv(SCREEN/"metrics/all_cases_metrics.csv")
    for row in metrics.itertuples(): render_case(SCREEN/"cases"/row.case_id)
    mosaics(); contact_sheet(metrics); gallery(metrics); review_template(metrics)


if __name__=="__main__":main()
