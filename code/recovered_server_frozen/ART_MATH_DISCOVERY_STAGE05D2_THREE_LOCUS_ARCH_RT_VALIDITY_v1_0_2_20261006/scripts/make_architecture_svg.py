#!/usr/bin/env python3
import sys, html
from pathlib import Path
import pandas as pd
p=Path(sys.argv[1]); out=Path(sys.argv[2])
a=pd.read_csv(p,sep='\t')
W,H=1400,760
m=80
rows=[]
for i,r in a.iterrows():
    lo=min(r.candidate_start0,r.repeat_envelope_start0,r.rt_start0)-600
    hi=max(r.candidate_end0,r.repeat_envelope_end0,r.rt_end0)+600
    y=150+i*190
    def x(v): return m+(v-lo)/(hi-lo)*(W-2*m)
    rows.append((r,y,x,lo,hi))
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
     '<rect width="100%" height="100%" fill="white"/>',
     '<style>text{font-family:Arial,sans-serif;fill:#111}.t{font-size:24px;font-weight:700}.s{font-size:17px}.n{font-size:14px}</style>',
     '<text x="70" y="50" class="t">Stage05D2: DNA-first three-locus architecture audit</text>',
     '<text x="70" y="80" class="s">Candidate = frozen DNA math region; repeat = prediction-relative noncoding repeat envelope; RT-like ORF = post-discovery protein interpretation</text>']
for r,y,x,lo,hi in rows:
    svg += [f'<text x="70" y="{y-55}" class="s">{html.escape(r.candidate_key)} · rank {int(r.region_rank)} · {html.escape(r.myrt_model)}</text>',
            f'<line x1="70" y1="{y}" x2="1330" y2="{y}" stroke="#777" stroke-width="2"/>',
            f'<rect x="{x(r.candidate_start0):.1f}" y="{y-24}" width="{x(r.candidate_end0)-x(r.candidate_start0):.1f}" height="48" fill="#d9eaf7" stroke="#1f5a7a"/>',
            f'<text x="{x(r.candidate_start0):.1f}" y="{y-31}" class="n">math candidate [{int(r.candidate_start0)},{int(r.candidate_end0)})</text>',
            f'<rect x="{x(r.repeat_envelope_start0):.1f}" y="{y+35}" width="{max(2,x(r.repeat_envelope_end0)-x(r.repeat_envelope_start0)):.1f}" height="24" fill="#ffd9b3" stroke="#a65b00"/>',
            f'<text x="{x(r.repeat_envelope_start0):.1f}" y="{y+80}" class="n">repeat envelope · gap {int(r.repeat_to_rt_gap_bp)} bp · {html.escape(r.repeat_side_relative_to_rt)}</text>',
            f'<line x1="{x(r.rt_start0):.1f}" y1="{y-42}" x2="{x(r.rt_end0):.1f}" y2="{y-42}" stroke="#722f8b" stroke-width="14"/>']
    # strand arrow
    if r.rt_strand=='-':
        svg += [f'<polygon points="{x(r.rt_start0):.1f},{y-42} {x(r.rt_start0)+18:.1f},{y-54} {x(r.rt_start0)+18:.1f},{y-30}" fill="#722f8b"/>']
    else:
        svg += [f'<polygon points="{x(r.rt_end0):.1f},{y-42} {x(r.rt_end0)-18:.1f},{y-54} {x(r.rt_end0)-18:.1f},{y-30}" fill="#722f8b"/>']
    svg += [f'<text x="{x(r.rt_start0):.1f}" y="{y-58}" class="n">RT-like ORF [{int(r.rt_start0)},{int(r.rt_end0)}) {html.escape(r.rt_strand)} · {int(r.rt_aa_len)} aa</text>',
            f'<text x="70" y="{y+112}" class="n">HMM model cov={r.myrt_model_coverage:.3f}; query cov={r.myrt_query_coverage:.3f}; YxDD={html.escape(str(r.yxdd_motifs))}</text>']
svg += ['<text x="70" y="720" class="n">Scope: biological interpretation only. No mathematical rescoring; no new-ART claim.</text>','</svg>']
out.write_text('\n'.join(svg))
