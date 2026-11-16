"""Create a local, self-contained HTML result report from real solver outputs."""
from pathlib import Path
import base64
from datetime import datetime
from zoneinfo import ZoneInfo
import html
import json


def build_report(folder):
    folder=Path(folder)
    sections=[]
    for case,title,image_name in [('prototype','Prototype 电力系统','prototype_dashboard'),
                                  ('softwarex','SoftwareX 工业转型','softwarex_dashboard')]:
        path=folder/case/'summary.json'
        if not path.exists(): continue
        s=json.loads(path.read_text(encoding='utf8'))
        if case=='prototype':
            content=f"最优解；目标值 {s['objective']:,.9f} 十亿美元；累计排放 {s['total_co2']:,.3f} Mt CO₂。"
            details=f"最大等式残差 {s['max_equality_residual']:.3g}；最大需求缺口 {s['max_demand_violation']:.3g} TWh。"
        else:
            content=f"{s['status']}；可行方案目标值 {s['objective_MUSD']:,.6f} 百万美元；最优性差距 {100*s['mip_gap']:.2f}%。"
            details=f"累计排放 {s['total_co2_tonnes']:,.6f} 吨 CO₂；全期预算 {s['cumulative_co2_budget_tonnes']:,.6f} 吨；最大约束残差 {s['max_constraint_violation']:.3g}（模型缩放单位）。"
        image=folder/case/'figures'/f'{image_name}.png'
        embed=''
        if image.exists():
            data=base64.b64encode(image.read_bytes()).decode('ascii')
            embed=f'<img alt="{title}结果图" src="data:image/png;base64,{data}">'
        sections.append(f'<section><h2>{title}</h2><p class="metric">{html.escape(content)}</p><p>{html.escape(details)}</p>{embed}<p><a href="{case}/summary.json">完整求解摘要</a> · <a href="{case}/solution_variables.csv">全部变量 CSV</a> · <a href="{case}/solver.log">求解日志</a></p></section>')
    timestamp=datetime.now(ZoneInfo('Asia/Shanghai')).strftime('%Y-%m-%d %H:%M %Z')
    document='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>STREAM Python 双案例结果</title><style>body{margin:0;background:#eef2f7;color:#253247;font:16px/1.7 "Microsoft YaHei",sans-serif}main{max-width:1250px;margin:auto;padding:32px}section{background:white;border-radius:16px;margin:24px 0;padding:24px;box-shadow:0 4px 20px #23334b0b}h1{font-size:32px}h2{font-size:23px}.metric{font-weight:bold}img{width:100%;height:auto}a{color:#356b99}small{color:#657489}</style><main><h1>STREAM · Python 重写与双案例求解</h1><p>Excel 输入 → Python 建模 → HiGHS → CSV 与 Matplotlib 图表</p>'''+''.join(sections)+f'''<section><h2>如何理解结果</h2><p>Prototype 取得最优解。SoftwareX 是混合整数问题；若显示 Time limit reached，图表展示经约束检查的可行方案，最优性差距仍存在，不能视为全局最优路线。</p><p>Prototype 使用 GW、TWh、Mt CO₂。SoftwareX 根据生成数据时的缩放定义，还原为吨、MMBtu 与百万美元。两个案例的排放单位不同。</p><p>图表与数字读取实际求解输出。SoftwareX 默认仅用原案例的 544 个离散决策作初始候选，没有固定这些变量。Julia 仅在迁移开发时用于参考系数核对，Python 日常运行无需 Julia。</p></section><small>生成时间：{timestamp}；报告图片已内嵌，可离线查看。</small></main></html>'''
    (folder/'report.html').write_text(document,encoding='utf8')
    return folder/'report.html'
