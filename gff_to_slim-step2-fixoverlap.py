#!/usr/bin/env python3
import argparse
import re

def parse_line(line):
    """
    从类似
      initializeGenomicElement(g1, 974898773, 974899147);
    的字符串中提取出 (feature, start, end) 三元组。
    """
    pattern = r'initializeGenomicElement\(\s*([^(,]+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\);'
    m = re.match(pattern, line)
    if m:
        feature = m.group(1).strip()
        start = int(m.group(2))
        end = int(m.group(3))
        return (feature, start, end)
    return None

def read_segments(filename):
    """
    读取输入文件（第一步输出结果），解析每一行返回区段列表，
    格式为 [(feature, start, end), ...]
    """
    segments = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            seg = parse_line(line)
            if seg:
                segments.append(seg)
    return segments

def fix_overlaps(segments):
    """
    第一轮修正：
    遍历已按起始位置排序的区段列表：
      - 若当前区段 cur 的起始位置 <= 前一区段 prev 的结束，则尝试将 cur.start
        调整为 prev.end + 1。
      - 如果调整后 cur.start < cur.end，保留调整结果；否则删除该区段。
    返回修正后的区段列表。
    """
    fixed = []
    segments_sorted = sorted(segments, key=lambda seg: seg[1])
    for seg in segments_sorted:
        if not fixed:
            fixed.append(seg)
        else:
            prev = fixed[-1]
            if seg[1] <= prev[2]:
                new_start = prev[2] + 1
                if new_start < seg[2]:
                    fixed.append((seg[0], new_start, seg[2]))
                # 否则跳过该区段
            else:
                fixed.append(seg)
    return fixed

def second_fix(segments):
    """
    第二轮修正：
    再次扫描第一轮修正后的列表，对仍存在重叠情况进行进一步调整：
      - 如果当前区段 cur 与前一区段 prev 重叠，则再次调整 cur.start = prev.end + 1，
        并检查是否有效；如无效则删除。
      - 此外可以添加针对“区段包裹”的判断（这里采用简单处理策略）。
    返回二次修正后的区段列表。
    """
    if not segments:
        return segments
    fixed = [segments[0]]
    for i in range(1, len(segments)):
        cur = segments[i]
        prev = fixed[-1]
        if cur[1] <= prev[2]:
            new_start = prev[2] + 1
            if new_start < cur[2]:
                fixed.append((cur[0], new_start, cur[2]))
            # 否则跳过 cur
        else:
            fixed.append(cur)
    return fixed

def output_segments(segments):
    """
    输出修正后的区段，格式为：
      initializeGenomicElement(feature, start, end);
    """
    for seg in segments:
        feature, start, end = seg
        print(f"initializeGenomicElement({feature}, {start}, {end});")

def main():
    parser = argparse.ArgumentParser(
        description="对第一步输出的 SLiM 区段结果进行二次修复（fix overlap）。"
    )
    parser.add_argument("-i", "--input", required=True, help="输入第一步输出的结果文件")
    args = parser.parse_args()
    
    # 读取输入区段
    segments = read_segments(args.input)
    if not segments:
        print("未能读取到有效的区段，请检查输入文件格式。")
        return

    # 第一轮修正
    fixed1 = fix_overlaps(segments)
    # 第二轮修正
    fixed2 = second_fix(fixed1)
    # 输出最终修正结果
    output_segments(fixed2)

if __name__ == "__main__":
    main()
