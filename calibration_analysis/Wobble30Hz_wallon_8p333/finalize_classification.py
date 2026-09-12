from pathlib import Path
HERE=Path(__file__).resolve().parent
p=Path(__file__).resolve().parent/'Wobble30Hz_wallon_8p333_rotational_load_report.md'
s=p.read_text(encoding='utf-8')
s=s.replace('## Decision\nB_SUSTAINED_POST_IMPACT_NONMAGNETIC_ROTATIONAL_STALL','## Decision\nE_POST_IMPACT_MECHANISM_STILL_AMBIGUOUS')
s += '\n\n## Reconstruction gate correction\nThe ODB does contain direct VR1-VR3. Relative-rotation-log and direct VR differ by RMS 51.392 rad/s and peak relative difference 0.992, exceeding the 10% consistency gate. Therefore the formal mechanism classification is conservatively E (ambiguous), despite the provisional window ratios suggesting nonmagnetic loading. The next action is to resolve the angular-velocity convention/coordinate basis before changing any physical parameter.\n'
p.write_text(s,encoding='utf-8')
(HERE/'mechanism_classification.csv').write_text('classification,reason\nE_POST_IMPACT_MECHANISM_STILL_AMBIGUOUS,Direct VR vs relative-rotation mismatch exceeds 10 percent\n',encoding='utf-8')
