#!/usr/bin/env python
# -*- coding:utf-8 -*-
import os
import re
import sys
import csv
import argparse
from collections import defaultdict
from jinja2 import (
    FileSystemLoader,
    Environment)
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8')

class Config(object):
    script_root_path = os.path.realpath(os.path.dirname(sys.argv[0]))
    static_path = os.path.join(script_root_path, "static")
    template_path = os.path.join(static_path, "template")


def initialized_environment():
    # 创建一个加载器, jinja2 会从这个目录中加载模板
    loader = FileSystemLoader(Config.template_path)
    # 用加载器创建一个环境, 有了它才能读取模板文件
    e = Environment(loader=loader)
    return e


def load_file(file):
    with open(file, 'r') as f:
        return f.read()


class Template(object):
    e = initialized_environment()

    @classmethod
    def render(cls, filename, *args, **kwargs):
        # 调用 get_template() 方法加载模板并返回
        template = cls.e.get_template(filename)
        # 用 render() 方法渲染模板
        # 可以传递参数
        return template.render(*args, **kwargs)


class TempFileNameBase(object):
    _keys = []

    def __init__(self, folder):
        self.folder = folder

    def readFile(self, file, callback=None):
        if callback is None:
            callback = self.default_callback
        return callback(file)

    def default_callback(self, file):
        r = {}
        for line in open(file):
            line = line.strip()
            k, v = line.split(',', 1)
            if v.isdigit() and "time" not in k:
                v = '{:,}'.format(int(v))
            k = k.replace(' ', '')
            r[k] = v
        return r


class TempFileNameSingle(TempFileNameBase):
    _keys = [
           "alignment_report",
           "annotated_report",
           "cell_report",
           "cluster",
           "cutoff",
           "marker",
           "sample",
           "sequencing_report",
           "vln"
    ]

    def getData(self):
        data = {}
        for key in self._keys:
            abspath = '{}.csv'.format(os.path.join(self.folder, key))
            print("开始读取文件[{}]".format(key))
            if key == 'vln':
                data[key] = self.readFile(abspath, self.vln_callback)
            elif key == 'cluster':
                maxNumi, minNumi, clusterMap, d = self.readFile(abspath, self.cluster_callback)
                data[key] = d
                data['maxNumi'] = maxNumi
                data['minNumi'] = minNumi
                data['clusterMap'] = clusterMap
            elif key == 'marker':
                data[key] = self.readFile(abspath, self.marker_callback)
            elif key == 'cutoff':
                data['cutoffMap'] = self.readFile(abspath, self.cutoff_callback)
            else:
                data[key] = self.readFile(abspath)
        return data

    def marker_callback(self, file):
        f = open(file)
        header = f.readline().strip()
        headers = [i.replace('.', '_') for i in header.split(',')]
        headers[0] = 'id'
        csv_dict = csv.DictReader(f, fieldnames=headers)
        d = []
        for item in csv_dict:
            item['p_val_adj'] = float(item['p_val_adj'])
            d.append(item)
        return d

    def vln_callback(self, file):
        d = {}
        f = open(file)
        header = f.readline().strip()
        headers = header.split(',')
        d[headers[1]] = []
        d[headers[2]] = []
        d[headers[3]] = []
        for line in f:
            name, a, b, c = line.strip().split(',')
            d[headers[1]].append(a)
            d[headers[2]].append(b)
            d[headers[3]].append(c)
        return d

    def cluster_callback(self, file):
        f = open(file)
        header = f.readline().strip()
        headers = header.split(',')
        headers[0] = 'id'
        csv_dict = csv.DictReader(f, fieldnames=headers)
        d = [dict(i) for i in csv_dict]
        maxNumi = sorted([int(float(item['nUMI'])) for item in d], reverse=True)[10]
        minNumi = min([int(float(item['nUMI'])) for item in d])
        clusterMap = defaultdict(list)
        for item in d:
            clusterMap[item['Cluster']].append(item)
        return maxNumi, minNumi, clusterMap, d

    def cutoff_callback(self, file):
        f = open(file)
        header = f.readline().strip()
        headers = header.split(',')
        csv_dict = csv.DictReader(f, fieldnames=headers)
        d = [dict(i) for i in csv_dict]
        cutoffMap = defaultdict(list)
        for item in d:
            cutoffMap[item['cell']].append(item)
        return cutoffMap


class TempFileNameTwo(TempFileNameBase):
    _keys = [
           "alignment_report",
           "annotated_report",
           "cell_report",
           "cutoff",
           "mix_report",
           "sample",
           "sequencing_report",
           "vln"
        ]

    def getStatic(self):
        path = Config.static_path
        js_files = []
        css_files = []
        for (root, dirs, files) in os.walk(path):
            for file in files:
                if file.endswith("css"):
                    css_files.append(os.path.join(root, file))
                elif file.endswith("js"):
                    js_files.append(os.path.join(root, file))
        js_string = '\n'.join((load_file(js_file) for js_file in js_files))
        css_string = '\n'.join((load_file(css_file) for css_file in css_files))
        return {"js": js_string, "css": css_string}

    def getData(self):
        data = {
            "static": self.getStatic()
        }
        for key in self._keys:
            abspath = '{}.csv'.format(os.path.join(self.folder, key))
            print("开始读取文件[{}]".format(key))
            if key == 'cutoff':
                data['cutoffMap'] = self.readFile(abspath, self.cutoff_callback)
            elif key == 'vln':
                data[key] = self.readFile(abspath, self.vln_callback)
            elif key == "annotated_report":
                data[key] = self.readFile(abspath, self.annotated_report_callback)
            elif key == 'cell_report':
                data['cell_report'] = self.readFile(abspath, self.cell_report_callback)
            else:
                data[key] = self.readFile(abspath)
        return data

    def vln_callback(self, file):
        d = {}
        f = open(file)
        header = f.readline().strip()
        headers = header.split(',')
        samples = []
        for v in headers:
            if v.endswith("UB") and v.split('_')[0] not in samples:
                samples.append(v.split('_')[0])
#        samples = list(set([v.split('_')[0] for v in headers if v.endswith("UB") and v.split('_')[0] not in ]))
        if samples.__len__() != 2:
            raise Exception("samples number != 2")
        sampleA = samples[0]
        sampleB = samples[1]
        d["Raw"] = [
            [], []
        ]
        d["UMI"] = [
            [], []
        ]
        d["Gene"] = [
            [], []
        ]

        d["Samples"] = samples

        d["mix"] = {
            "sampleA": [],
            "sampleB": [],
            "mix": []
        }
        for line in f:
            name, Raw, SampleA_UB, SampleA_GN, SampleB_UB, SampleB_GN, Species_UMI = line.strip().split(',')
            if Species_UMI != "Mix":
                if Species_UMI == sampleA:
                    d["Raw"][0].append(Raw)
                    d["UMI"][0].append(SampleA_UB)
                    d["Gene"][0].append(SampleA_GN)
                    d["mix"]["sampleA"].append([SampleA_UB, SampleB_UB])
                elif Species_UMI == sampleB:
                    d["Raw"][1].append(Raw)
                    d["UMI"][1].append(SampleB_UB)
                    d["Gene"][1].append(SampleB_GN)
                    d["mix"]["sampleB"].append([SampleA_UB, SampleB_UB])
            else:
                d["mix"]["mix"].append([SampleA_UB, SampleB_UB])
        return d

    def cutoff_callback(self, file):
        f = open(file)
        header = f.readline().strip()
        headers = header.split(',')
        csv_dict = csv.DictReader(f, fieldnames=headers)
        d = [dict(i) for i in csv_dict]
        cutoffMap = defaultdict(list)
        for item in d:
            cutoffMap[item['cell']].append(item)
        return cutoffMap

    def cell_report_callback(self, file):
        file_string = load_file(file)
        samples = list(set(re.findall(r"\((.+?)\)", file_string)))
        if samples.__len__() != 2:
            raise Exception("samples number != 2")
        # 原始的map
        data_map = self.default_callback(file)
        r = {}
        sampleA = samples[0]
        sampleB = samples[1]
        r["sampleA"] = {"name": sampleA}
        r["sampleB"] = {"name": sampleB}
        for k, v in data_map.items():
            if sampleA in k:
                k = k.replace('(', '').replace(')', '').replace(sampleA, '')
                r["sampleA"][k] = v
            elif sampleB in k:
                k = k.replace('(', '').replace(')', '').replace(sampleB, '')
                r["sampleB"][k] = v
            else:
                r[k] = v
        return r

    def annotated_report_callback(self, file):
        file_string = load_file(file)
        samples = list(set(re.findall(r"\((.+?)\)", file_string)))
        data_map = self.default_callback(file)
        r = {}
        if samples.__len__() == 0:
            sampleA = 'GRCh38'
            sampleB = 'mm10'
            r={'sampleA': {'name': 'GRCh38', 'ReadsMappedConfidentlytoGenome': '-', 'ReadsMappedConfidentlytoGene': '-', 'ReadsMappedConfidentlytoExonicRegions': '-', 'ReadsMappedConfidentlytoIntronicRegions':     '-', 'ReadsMappedAntisensetoGene': '-'}, 'sampleB': {'name': 'mm10', 'ReadsMappedConfidentlytoGenome': '-', 'ReadsMappedConfidentlytoGene': '-', 'ReadsMappedConfidentlytoExonicRegions': '-', 'ReadsMapped    ConfidentlytoIntronicRegions': '-', 'ReadsMappedAntisensetoGene': '-'}}
        elif samples.__len__() == 2:
            sampleA = samples[0]
            sampleB = samples[1]
            r["sampleA"] = {"name": sampleA}
            r["sampleB"] = {"name": sampleB}
        else:
            raise Exception("samples number != 2")
        # 原始的map
        for k, v in data_map.items():
            if sampleA in k:
                k = k.replace('(', '').replace(')', '').replace(sampleA, '')
                r["sampleA"][k] = v
            elif sampleB in k:
                k = k.replace('(', '').replace(')', '').replace(sampleB, '')
                r["sampleB"][k] = v
            else:
                r[k] = v
        return r


def get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='commands')
    # single
    sub_single = subparsers.add_parser("single", help="单样本")
    sub_single.add_argument(dest='folder', help="输入目录")
    sub_single.add_argument(dest='out_prefix', help="输出文件前缀")
    # 绑定函数
    sub_single.set_defaults(func=single)
    # two
    sub_two = subparsers.add_parser("two", help="双样本")
    sub_two.add_argument(dest='folder', help="输入目录")
    sub_two.add_argument(dest='out_prefix', help="输出文件前缀")
    # 绑定函数
    sub_two.set_defaults(func=two)
    # 如果没有输入-h参数则在命令行参数上添加-h
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    return parser.parse_args()


# single 入口函数
def single(args):
    tf = TempFileNameSingle(args.folder)
    data = tf.getData()
    html = args.out_prefix+'.html'
    print('开始生成HTML文件[{}]'.format(html))
    path = os.path.join(args.folder, html)
    with open(path, 'w', encoding="utf-8") as f:
        print(Template.render('htmlTemp.single.html', **data, n=''), file=f)


# two 入口函数
def two(args):
    tf = TempFileNameTwo(args.folder)
    data = tf.getData()
    html = args.out_prefix+'.html'
    print('开始生成HTML文件[{}]'.format(html))
    path = os.path.join(args.folder, html)
    with open(path, 'w', encoding="utf-8") as f:
        print(Template.render('htmlTemp.two.html', **data, n=''), file=f)


# 主程序
def main():
    # 获取参数, 子命令绑定入口函数
    args = get_args()
    # 得到参数对象并调用绑定的函数, 传入参数
    # single(args) or two(args)
    args.func(args)


if __name__ == '__main__':
    main()
