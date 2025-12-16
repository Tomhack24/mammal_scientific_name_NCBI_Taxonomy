#!/usr/bin/env python3
"""
NCBI Taxonomy データベースから哺乳類（Mammalia）の確定した学名のみを取得するスクリプト
（sp., cf., aff., x（ハイブリッド）、その他の暫定的な名前を除外）
"""

import os
import re
from collections import defaultdict

# 設定
TAXDUMP_DIR = "taxdump"
MAMMALIA_TAX_ID = 40674


def is_valid_species_name(name):
    """
    確定した学名かどうかを判定
    以下を除外:
    - sp. を含む（未同定種）
    - cf. を含む（比較参照）
    - aff. を含む（近縁）
    - x を含む（ハイブリッド）
    - 'isolate'、'haplotype'などの文字列を含む
    - カッコを含む（provisional names）
    - BOLDなどの識別子を含む
    """
    name_lower = name.lower()
    
    # 除外パターン
    exclude_patterns = [
        r'\bsp\.',           # sp. (species)
        r'\bcf\.',           # cf. (confer)
        r'\baff\.',          # aff. (affinis)
        r'\bx\b',            # x (hybrid marker)
        r'\bisolate\b',      # isolate
        r'\bhaplotype\b',    # haplotype
        r'\bgen\.',          # gen. (genus)
        r'\bgen\b',          # gen
        r'\benvironmental\b', # environmental sample
        r'\bmixed\b',        # mixed
        r'\bindet\.',        # indeterminate
        r'bold:',            # BOLD identifier
        r'[\(\)]',           # parentheses
        r'\d{4}',            # 4-digit numbers (often year codes)
        r'[A-Z]{2,}-\d{4}',  # codes like AG-2015
        r'[A-Z]{3,}\d+',     # codes like FMNH129874
        r'\b(complex|group)\b', # complex or group
    ]
    
    for pattern in exclude_patterns:
        if re.search(pattern, name_lower):
            return False
    
    # 学名は通常2語（属名 + 種小名）または3語（属名 + 種小名 + 亜種名）
    # それ以上の単語がある場合は暫定的な名前の可能性が高い
    words = name.split()
    
    # 2語または3語の学名のみを受け入れる
    if len(words) < 2 or len(words) > 3:
        return False
    
    # 最初の文字は大文字、2番目以降は小文字で始まるべき
    if not words[0][0].isupper():
        return False
    
    for word in words[1:]:
        if not word[0].islower():
            return False
    
    return True


def load_nodes(nodes_file):
    """nodes.dmpを読み込み、親子関係とrankを辞書に格納"""
    print(f"Loading {nodes_file}...")
    parent_dict = {}  # tax_id -> parent_tax_id
    rank_dict = {}    # tax_id -> rank
    children_dict = defaultdict(list)  # parent_tax_id -> [child_tax_ids]
    
    with open(nodes_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t|\t')
            tax_id = int(parts[0].strip())
            parent_tax_id = int(parts[1].strip())
            rank = parts[2].strip()
            
            parent_dict[tax_id] = parent_tax_id
            rank_dict[tax_id] = rank
            children_dict[parent_tax_id].append(tax_id)
    
    print(f"Loaded {len(parent_dict)} nodes")
    return parent_dict, rank_dict, children_dict


def load_names(names_file):
    """names.dmpを読み込み、scientific nameを取得"""
    print(f"Loading {names_file}...")
    names_dict = {}  # tax_id -> scientific name
    
    with open(names_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t|\t')
            tax_id = int(parts[0].strip())
            name = parts[1].strip()
            name_class = parts[3].strip().rstrip('\t|')
            
            # scientific nameのみ取得
            if name_class == "scientific name":
                names_dict[tax_id] = name
    
    print(f"Loaded {len(names_dict)} scientific names")
    return names_dict


def get_all_descendants(tax_id, children_dict):
    """指定されたtax_idの全ての子孫を再帰的に取得"""
    descendants = set()
    queue = [tax_id]
    
    while queue:
        current = queue.pop(0)
        if current in children_dict:
            for child in children_dict[current]:
                if child not in descendants:
                    descendants.add(child)
                    queue.append(child)
    
    return descendants


def main():
    # データファイルのパス
    nodes_file = os.path.join(TAXDUMP_DIR, "nodes.dmp")
    names_file = os.path.join(TAXDUMP_DIR, "names.dmp")
    
    # データ読み込み
    parent_dict, rank_dict, children_dict = load_nodes(nodes_file)
    names_dict = load_names(names_file)
    
    # Mammaliaの全子孫を取得
    print(f"\nFinding all descendants of Mammalia (tax_id: {MAMMALIA_TAX_ID})...")
    descendants = get_all_descendants(MAMMALIA_TAX_ID, children_dict)
    print(f"Found {len(descendants)} descendants")
    
    # 種レベル（rank = "species"）のものを抽出
    print("\nFiltering species...")
    all_mammal_species = []
    valid_mammal_species = []
    
    for tax_id in descendants:
        if rank_dict.get(tax_id) == "species":
            name = names_dict.get(tax_id, "Unknown")
            all_mammal_species.append(name)
            
            # 確定した学名のみを抽出
            if is_valid_species_name(name):
                valid_mammal_species.append(name)
    
    # 学名でソート
    all_mammal_species.sort()
    valid_mammal_species.sort()
    
    # 結果を表示
    print(f"\nTotal mammal species: {len(all_mammal_species)}")
    print(f"Valid species names (confirmed): {len(valid_mammal_species)}")
    print(f"Excluded (provisional names): {len(all_mammal_species) - len(valid_mammal_species)}")
    
    print("\nFirst 20 valid mammal species:")
    for i, name in enumerate(valid_mammal_species[:20], 1):
        print(f"{i:4d}. {name}")
    
    # 結果をファイルに保存（確定した学名のみ）
    output_file = "mammal_scientific-nameb.txt"
    with open(output_file, 'w', encoding='utf-8') as f:
        for name in valid_mammal_species:
            f.write(f"{name}\n")
    
    print(f"\nResults saved to {output_file}")
    print(f"Total confirmed species names: {len(valid_mammal_species)}")


if __name__ == "__main__":
    main()
