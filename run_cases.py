"""Command line: python run_cases.py --case all (no Julia required)."""
from pathlib import Path
import argparse
import json


def main():
    root=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description='STREAM Python / HiGHS case studies')
    parser.add_argument('--case',choices=['all','prototype','softwarex'],default='all')
    parser.add_argument('--output',type=Path,default=root/'results')
    parser.add_argument('--prototype-input',type=Path,default=root/'inputs/prototype/prototype_0.xlsx')
    parser.add_argument('--softwarex-input',type=Path,default=root/'inputs/softwarex/f_4.xlsx')
    parser.add_argument('--time-limit',type=float,default=300)
    parser.add_argument('--mip-gap',type=float,default=.05)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--warm-start-dir',type=Path,default=root/'inputs/softwarex/warmstart')
    parser.add_argument('--cold-start',action='store_true',help='Ignore saved binary incumbent')
    parser.add_argument('--plots-only',action='store_true')
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args()
    if not args.plots_only:
        if args.case in ['all','prototype']:
            from stream_python.prototype import run_prototype
            r=run_prototype(args.prototype_input,args.output/'prototype')
            print(json.dumps(r,indent=2,ensure_ascii=False))
        if args.case in ['all','softwarex']:
            from stream_python.softwarex import run_softwarex
            r=run_softwarex(args.softwarex_input,args.output/'softwarex',
                            time_limit=args.time_limit,mip_gap=args.mip_gap,threads=args.threads,
                            warm_start_dir=None if args.cold_start else args.warm_start_dir)
            print(json.dumps(r,indent=2,ensure_ascii=False))
    if not args.no_plots:
        from stream_python.plots import plot_all
        plot_all(args.output)
        from stream_python.report import build_report
        build_report(args.output)
    print(f'Results: {args.output.resolve()}')


if __name__=='__main__':
    main()
