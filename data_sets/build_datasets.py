import pandas as pd, glob, re
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment, PatternFill
EX='../meta_analysis/extractions/'
# study lookup
rows=[]
for f in glob.glob(EX+'*.csv'):
    try: d=pd.read_csv(f)
    except Exception: continue
    if 'study_id' in d.columns and 'doi' in d.columns:
        for c in ['authors','year','title','journal','doi']:
            if c not in d.columns: d[c]=''
        rows.append(d[['study_id','authors','year','title','journal','doi']])
L=pd.concat(rows).dropna(subset=['study_id']).drop_duplicates('study_id').set_index('study_id')
T=pd.read_csv('../meta_analysis/synthesis/meta_table.csv')
def lk(s,c):
    return L.loc[s,c] if s in L.index else ''
for c in ['authors','year','title','journal','doi']: T[c]=T.study.map(lambda s: lk(s,c))
T.loc[T.study.str.startswith('GLORICH'),['authors','year','title','journal','doi']]=['This study, from Hartmann et al. (2014) GLORICH','','Within-station mixed-effects model (Text S2)','','']
T.loc[T.cls.str.contains('HFO|NICA'),['authors','year','title','journal','doi']]=['This study','','Speciation model run (Materials and methods)','','']
pub=~T.cls.str.contains('HFO|NICA')&~T.study.str.startswith('GLORICH')
T.insert(0,'source',['Published' if p else 'This study' for p in pub])
T=T.rename(columns={'cls':'evidence_class','slope':'slope_dlog10C_dpH','se':'standard_error','lo':'lower_95','hi':'upper_95','note':'quality_or_note'})
P=pd.read_csv('../meta_analysis/synthesis/meta_pooled.csv')
Vt=pd.read_csv('../meta_analysis/synthesis/diel_votes.csv')
missing=sorted(set(T[pub & (T.authors.fillna('')=='')].study)); print('studies without citation:',missing)
readme1=pd.DataFrame({'Data Set S1':[
 'Hartland et al. Inland waters are not small oceans. Supporting Information, Data Set S1.',
 'Published responses of labile, free and dissolved trace metals and oxyanions to pH in fresh water, with speciation model results.',
 f'Sheet "Slopes": {pub.sum()} published slopes from {T[pub].study.nunique()} studies, plus slopes from this study (GLORICH) and speciation model results.',
 'slope_dlog10C_dpH: change in log10 concentration per pH unit. standard_error, lower_95, upper_95: as reported or computed from tabulated data (blank where not available).',
 'evidence_class: how well the design isolates pH (see Materials and methods and Fig. 3c). The electrostatic limit is -2 for divalent cations; phi = -slope/2.',
 'Sheet "Summary": median, interquartile range and range by element and group, and by evidence class for divalent cations other than Fe; phi columns give the fraction of the electrostatic limit.',
 'Sheet "Diel votes": number of diel studies in which each element fell, did not change or rose by day.',
 'Omanovic et al. (2015) slopes reflect mixing of mine water and are listed but not summarised.']})
with pd.ExcelWriter('Data_Set_S1_meta_analysis.xlsx') as w:
    readme1.to_excel(w,sheet_name='README',index=False); T.to_excel(w,sheet_name='Slopes',index=False); P.to_excel(w,sheet_name='Summary',index=False); Vt.to_excel(w,sheet_name='Diel votes',index=False)
# ---- S2
def rd(f): return pd.read_csv(EX+f)
M=pd.read_csv('../meta_analysis/synthesis/metab_rows.csv')
GE=pd.concat([rd('metab_gpp_er.csv'),rd('metab2_gpp_er.csv')],ignore_index=True)
sheets={'Metabolism Ea':M,'Metabolism GPP-ER':GE,'Oxygen and redox':rd('redox3_rows.csv'),'Carbon uptake':rd('ccm3_rows.csv'),'Food webs':rd('food3_rows.csv'),'Invaders':rd('invasive3_rows.csv')}
ox=pd.DataFrame([
 dict(study_id='Nkwalale2026',authors='Nkwalale LGT et al.',year=2026,journal='Communications Sustainability 1:86',doi='10.1038/s44458-026-00093-z',variable='Lake anoxia projections, 73 lakes, 2015-2099',value='SSP5-8.5: anoxia onset about 30 days earlier; anoxic share in eutrophic lakes about 40% to 60% by 2099',location_in_paper='pp. 1-2'),
 dict(study_id='Guan2026',authors='Guan Q, Shi K, Pi X',year=2026,journal='Science Advances 12:eaef3132',doi='10.1126/sciadv.aef3132',variable='River dissolved oxygen trend, 21,439 reaches, 1985-2023',value='-0.045 mg/L per decade; 78.8% of reaches deoxygenating; 21st-century decline 1.1% (SSP1-2.6) to 4.7% (SSP5-8.5)',location_in_paper='pp. 1, 4'),
 dict(study_id='Wang2025',authors='Wang J et al.',year=2025,journal='Science Advances 11:eadr1695',doi='10.1126/sciadv.adr1695',variable='Global inland-water oxygen budget 1900-2010',value='Production 0.16 to 0.94 Pg/yr; consumption 0.44 to 1.47 Pg/yr; warming raised consumption (+0.2 Pg/yr) more than production (+0.1 Pg/yr)',location_in_paper='pp. 1, 6'),
 dict(study_id='Zhang2026',authors='Zhang W, Zhang D, Shan B, Zhou Q',year=2026,journal='ACS ES&T Water 6:541-544',doi='10.1021/acsestwater.5c01231',variable='Share of lake deoxygenation from solubility',value='55% (cited)',location_in_paper='p. 1')])
sheets['Oxygen and redox']=pd.concat([sheets['Oxygen and redox'],ox],ignore_index=True)
readme2=pd.DataFrame({'Data Set S2':[
 'Hartland et al. Inland waters are not small oceans. Supporting Information, Data Set S2.',
 'Literature syntheses used in the framework and the Discussion. Each row gives the source, the value as reported and its location in the paper.',
 'Metabolism Ea: apparent activation energies of gross primary production and ecosystem respiration from within-system studies (Table S3).',
 'Metabolism GPP-ER: ratios of gross primary production to ecosystem respiration.',
 'Oxygen and redox: oxygen thresholds, sediment oxygen demand, release of metals and phosphorus at low oxygen, and global oxygen trends.',
 'Carbon uptake: half-saturation constants and free CO2 minima relevant to carbon-concentrating mechanisms.',
 'Food webs: hypoxia thresholds for New Zealand fauna and carbonate thresholds for shell builders.',
 'Invaders: effects of invasive bivalves and carp on oxygen, CO2, pH and nutrients.',
 'Values marked as computed were derived from published data; quality A, B or C reflects how directly the value is reported (A: tabulated or stated; C: read from a figure or derived).']})
with pd.ExcelWriter('Data_Set_S2_literature_syntheses.xlsx') as w:
    readme2.to_excel(w,sheet_name='README',index=False)
    for k,v in sheets.items(): v.to_excel(w,sheet_name=k[:31],index=False)
# styling
for f in ['Data_Set_S1_meta_analysis.xlsx','Data_Set_S2_literature_syntheses.xlsx']:
    p=''+f; wb=load_workbook(p)
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for c in row: c.font=Font(name='Arial',size=10,bold=(c.row==1))
        for c in ws[1]: c.fill=PatternFill('solid',fgColor='E8EEF4')
        ws.freeze_panes='A2'
        for col in ws.columns:
            w=max(len(str(c.value)) if c.value is not None else 0 for c in col[:200])
            ws.column_dimensions[col[0].column_letter].width=min(max(10,w+2),60 if ws.title!='README' else 140)
    wb.save(p)
for k,v in sheets.items(): print(k, v.shape)
print('S1',T.shape,pub.sum(),T[pub].study.nunique())
