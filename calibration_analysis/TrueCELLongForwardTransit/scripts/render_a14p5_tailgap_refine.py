"""Render fixed-view and matched-phase outputs for the TAIL-gap mesh test."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
sys.path.insert(0, str(HERE))
import analyze_a14p5_tailgap_refine as analysis
sys.path.insert(0, str(REPO / "calibration_analysis" / "FastPrecomputedForwardFlowScreen" / "scripts"))
import analyze_render_forward_flow as legacy


BLUE, RED, BODY = "#2864a8", "#d1493f", "#383c42"


def visual_fields(data):
    source = np.load(data["case"] / "private" / "truecel_field_private.npz")
    time = source["time"].astype(float)
    centers = source["fluid_element_centroids_mm"].astype(float)
    conn = source["fluid_connectivity_index"].astype(int)
    evf = np.clip(source["fluid_evf"].astype(float), 0.0, 1.0)
    pressure = source["fluid_pressure"].astype(float)
    velocity = np.nanmean(source["fluid_velocity_mm_s"][:, conn, :], axis=2)
    relative = centers - data["pipe"]
    return {"time":time, "centers":centers, "evf":evf, "pressure":pressure, "velocity":velocity,
            "s":relative@data["c"], "n":relative@data["n"], "b":relative@data["b"],
            "vs":velocity@data["c"], "vn":velocity@data["n"], "speed":np.linalg.norm(velocity,axis=2)}


def robot_geometry(data):
    deck=(data["case"]/f"{data['job']}.inp").read_text(encoding="latin1")
    _,nodes=legacy.part_nodes(deck,"Robot_SOLID")
    material_s=(nodes-data["center0"])@data["c"]
    colors=np.full(len(nodes),BODY,object)
    colors[material_s<=material_s.min()+.25]=BLUE
    colors[material_s>=material_s.max()-.25]=RED
    return nodes,colors


def robot_at(data,nodes,time):
    u=np.asarray([np.interp(time,data["time"],data["u"][:,j]) for j in range(3)])
    ur=np.asarray([np.interp(time,data["time"],data["ur"][:,j]) for j in range(3)])
    return data["center0"]+u+Rotation.from_rotvec(ur).apply(nodes-data["center0"])


def draw_robot(ax,data,nodes,colors,time,labels=True):
    points=robot_at(data,nodes,time)
    relative=points-data["pipe"]
    x,y=relative@data["c"],relative@data["n"]
    ax.scatter(x[::2],y[::2],c=colors[::2],s=7,linewidths=0,zorder=5)
    tail,head=int(np.argmin(x)),int(np.argmax(x))
    ax.scatter([x[tail],x[head]],[y[tail],y[head]],c=[BLUE,RED],s=30,zorder=6)
    if labels:
        ax.text(x[tail],y[tail]+.075,"TAIL",color=BLUE,ha="center",fontsize=7,weight="bold")
        ax.text(x[head],y[head]+.075,"HEAD",color=RED,ha="center",fontsize=7,weight="bold")
    return points[tail]


def side_slice(fields,k,b_value=0.0,width=.06):
    wet=fields["evf"][k]>.01
    mask=wet&(np.abs(fields["b"]-b_value)<=width)
    return mask


def fluid_scatter(ax,fields,k,quantity,vmin,vmax,b_value=0.0,zoom=False,vectors=False):
    mask=side_slice(fields,k,b_value,max(.055,width_for(fields,b_value)))
    values=quantity[k,mask]
    size=18 if zoom else 9
    image=ax.scatter(fields["s"][mask],fields["n"][mask],c=values,s=size,marker="s",
                     cmap="coolwarm",vmin=vmin,vmax=vmax,linewidths=0,rasterized=True)
    if vectors:
        ids=np.flatnonzero(mask)[::max(1,int(np.count_nonzero(mask)/80))]
        ax.quiver(fields["s"][ids],fields["n"][ids],fields["vs"][k,ids],fields["vn"][k,ids],
                  color="black",alpha=.65,angles="xy",scale_units="xy",scale=2500,width=.0022)
    return image


def width_for(fields,b_value):
    distances=np.unique(np.round(np.abs(fields["b"]-b_value),6))
    return float(distances[min(1,len(distances)-1)]+1e-6)


def common_grid(case,fields,time,b_value,quantity):
    k=int(np.argmin(np.abs(fields["time"]-time)))
    mask=side_slice(fields,k,b_value,max(.065,width_for(fields,b_value)))
    gx=np.linspace(-5.3,-3.1,150);gy=np.linspace(-.75,.75,100)
    xx,yy=np.meshgrid(gx,gy)
    values=griddata(np.column_stack((fields["s"][mask],fields["n"][mask])),quantity[k,mask],(xx,yy),
                    method="linear")
    return gx,gy,values


def main():
    coarse=analysis.load_case(analysis.COARSE);fine=analysis.load_case(analysis.FINE)
    coarse_fields=visual_fields(coarse);fine_fields=visual_fields(fine)
    coarse_nodes,coarse_colors=robot_geometry(coarse);fine_nodes,fine_colors=robot_geometry(fine)
    radius=float(fine["identity"]["lumen_radius_mm"])

    # Main fixed side-view robot GIF.
    frame_ids=np.unique(np.linspace(0,len(fine_fields["time"])-1,134).astype(int))
    fig,ax=plt.subplots(figsize=(11,3.7))
    def draw_main(frame):
        k=frame_ids[frame];t=fine_fields["time"][k];ax.clear()
        ax.axhline(radius,color="#5996a5");ax.axhline(-radius,color="#5996a5")
        draw_robot(ax,fine,fine_nodes,fine_colors,t)
        ds=np.interp(t,fine["time"],fine["axial"]-fine["axial"][0])
        vs=np.interp(t,fine["time"],fine["vs"])
        ax.set(xlim=(-6,6),ylim=(-.8,.8),aspect="equal",xlabel="canonical +s (mm), left to right",ylabel="n_routeA (mm)",
               title=f"{analysis.FINE} | t={t*1e3:.3f} ms | 120 Hz | B0=12 mT | G=2.0 mT\n"
                     f"delta_s={ds:+.5f} mm | v_s={vs:+.2f} mm/s | local radial r-refinement")
        ax.grid(alpha=.15);fig.tight_layout()
    FuncAnimation(fig,draw_main,frames=len(frame_ids),interval=80).save(
        OUT/f"{analysis.FINE}.gif",writer=PillowWriter(fps=12.5),dpi=100)
    plt.close(fig)

    # Whole fluid side slice with fixed pressure scale.
    finite=np.abs(fine_fields["pressure"][fine_fields["evf"]>.01])
    pclip=max(float(np.percentile(finite,99)),1e-9)
    fluid_ids=np.unique(np.linspace(0,len(fine_fields["time"])-1,120).astype(int))
    fig,ax=plt.subplots(figsize=(11,3.8))
    def draw_fluid(frame):
        k=fluid_ids[frame];t=fine_fields["time"][k];ax.clear()
        fluid_scatter(ax,fine_fields,k,fine_fields["pressure"],-pclip,pclip,0.0,False,True)
        ax.axhline(radius,color="#25282c");ax.axhline(-radius,color="#25282c")
        draw_robot(ax,fine,fine_nodes,fine_colors,t)
        ax.set(xlim=(-6,6),ylim=(-.8,.8),aspect="equal",xlabel="canonical +s (mm)",ylabel="n_routeA (mm)",
               title=f"{analysis.FINE} fluid | t={t*1e3:.3f} ms | pressure fixed scale +/-{pclip:.3g} N/mm2")
        fig.tight_layout()
    FuncAnimation(fig,draw_fluid,frames=len(fluid_ids),interval=80).save(
        OUT/"TRUECEL_A14P5_TAILGAP_REFINE_FLUID.gif",writer=PillowWriter(fps=12.5),dpi=100)
    plt.close(fig)

    # Dynamic TAIL-gap plane and zoom during the critical window.
    zoom_ids=np.flatnonzero((fine_fields["time"]>=.008)&(fine_fields["time"]<=.0125))[::2]
    fig,ax=plt.subplots(figsize=(10,4.2))
    def draw_zoom(frame):
        k=zoom_ids[frame];t=fine_fields["time"][k];ax.clear()
        tail=draw_robot(ax,fine,fine_nodes,fine_colors,t)
        tail_rel=tail-fine["pipe"]
        tail_s=float(tail_rel@fine["c"]);tail_b=float(tail_rel@fine["b"])
        fluid_scatter(ax,fine_fields,k,fine_fields["pressure"],-pclip,pclip,tail_b,True,True)
        draw_robot(ax,fine,fine_nodes,fine_colors,t)
        ax.axhline(radius,color="#25282c");ax.axhline(-radius,color="#25282c")
        ax.set(xlim=(tail_s-.75,tail_s+.75),ylim=(-.75,.75),aspect="equal",xlabel="canonical +s (mm)",ylabel="n_routeA (mm)",
               title=f"TAIL-gap flow zoom | t={t*1e3:.3f} ms | slice b={tail_b:+.3f} mm\npressure + velocity vectors; fixed pressure scale")
        fig.tight_layout()
    FuncAnimation(fig,draw_zoom,frames=len(zoom_ids),interval=90).save(
        OUT/"TRUECEL_A14P5_TAILGAP_FLOW_ZOOM.gif",writer=PillowWriter(fps=11.1),dpi=110)
    plt.close(fig)

    # Matched absolute times are also matched magnetic phases because f/start are frozen.
    closure=analysis.force_closure(fine)
    critical=(closure["history_time"]>=.008)&(closure["history_time"]<=.0125)
    peak=float(closure["history_time"][np.flatnonzero(critical)[np.argmin(closure["ffluid"][critical])]])
    tail_windows=fine["contact"]["TAIL"].get("episode_windows_ms",[])
    recontact=next((row[0]*1e-3 for row in tail_windows if row[0]>=analysis.T*1e3),.009275)
    phase_rows=(("pre-squeeze",max(.008,peak-.0001)),("peak negative fluid force",peak),
                ("post-peak",min(.0125,peak+.0001)),("TAIL recontact",recontact))
    variables=(("pressure",lambda f:f["pressure"],"N/mm2"),
               ("axial velocity",lambda f:f["vs"],"mm/s"),("EVF",lambda f:f["evf"],"-"))
    fig,axes=plt.subplots(len(phase_rows)*len(variables),3,figsize=(13,29),squeeze=False)
    row=0
    for phase_name,t in phase_rows:
        fine_tail=robot_at(fine,fine_nodes,t)[int(np.argmin((robot_at(fine,fine_nodes,t)-fine["pipe"])@fine["c"]))]
        b_value=float((fine_tail-fine["pipe"])@fine["b"])
        for var_name,getter,units in variables:
            gx,gy,a=common_grid(coarse,coarse_fields,t,b_value,getter(coarse_fields))
            _,_,z=common_grid(fine,fine_fields,t,b_value,getter(fine_fields))
            d=z-a
            finite_values=np.abs(np.r_[a[np.isfinite(a)],z[np.isfinite(z)]])
            vmax=max(float(np.percentile(finite_values,99)),1e-12)
            if var_name=="EVF": vmin,vmax_main,cmap=0,1,"viridis"
            else: vmin,vmax_main,cmap=-vmax,vmax,"coolwarm"
            diffmax=max(float(np.nanpercentile(np.abs(d),99)),1e-12)
            for col,(values,title,lo,hi,cm) in enumerate(((a,"original",vmin,vmax_main,cmap),
                                                         (z,"TAIL-refined",vmin,vmax_main,cmap),
                                                         (d,"difference",-diffmax,diffmax,"coolwarm"))):
                im=axes[row,col].imshow(values,origin="lower",extent=(gx[0],gx[-1],gy[0],gy[-1]),
                                        aspect="auto",cmap=cm,vmin=lo,vmax=hi)
                axes[row,col].axhline(radius,color="black",lw=.6);axes[row,col].axhline(-radius,color="black",lw=.6)
                axes[row,col].set_title(f"{phase_name} {t*1e3:.3f} ms | {var_name} | {title}",fontsize=8)
                axes[row,col].set_ylabel("n_routeA (mm)");fig.colorbar(im,ax=axes[row,col],fraction=.025,pad=.01,label=units)
            row+=1
    for ax in axes[-1]:ax.set_xlabel("canonical s (mm)")
    fig.suptitle("A14P5 coarse vs TAIL-gap local mesh at matched magnetic phases",y=.998)
    fig.tight_layout()
    fig.savefig(OUT/"TRUECEL_TAILREFINE_COARSE_VS_FINE_PHASE_COMPARE.png",dpi=170)
    plt.close(fig)


if __name__=="__main__":
    main()
