#!/usr/bin/env python3
"""Generate the public index and validate paper metadata (Python standard library)."""
import argparse
import hashlib
import json
import re
import tempfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


SLUG = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*")
HASH = re.compile(r"[0-9a-f]{64}")
WRAPPERS = {'paper.json', 'README.md', 'CHANGELOG.md', 'citation.bib'}
GENERATED_SUFFIXES = {'.pyc', '.aux', '.log', '.out', '.bbl', '.blg', '.fls', '.fdb_latexmk'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def text(value):
    return isinstance(value, str) and bool(value.strip())


def safe_file(base, relative):
    require(text(relative), 'artifact path must be a nonempty string')
    rel = Path(relative)
    require(not rel.is_absolute() and '..' not in rel.parts and
            rel.as_posix() == relative and '\\' not in relative,
            f'unsafe artifact path: {relative}')
    target = base / rel
    require(target.resolve().is_relative_to(base.resolve()),
            f'unsafe artifact path: {relative}')
    require(not base.is_symlink() and all(not (base / Path(*rel.parts[:i])).is_symlink()
            for i in range(1, len(rel.parts) + 1)), f'symlink artifact path: {relative}')
    require(target.is_file(), f'missing required file: {relative}')
    return target


def validate_paper(path, refresh=False):
    require(not path.is_symlink(), f'symlink metadata: {path}')
    paper = json.loads(path.read_text(encoding='utf-8'))
    require(isinstance(paper, dict), f'{path}: metadata must be an object')
    require(type(paper.get('schema_version')) is int and paper['schema_version'] == 1,
            f'{path}: schema_version must be 1')
    for key in ('slug', 'title', 'version', 'abstract'):
        require(text(paper.get(key)), f'{path}: {key} must be a nonempty string')
    require(SLUG.fullmatch(paper['slug']) is not None and paper['slug'] == path.parent.name,
            f'{path}: slug must match its lowercase, hyphenated directory name')
    require(re.fullmatch(r'\d{4}-\d{2}-\d{2}', paper['version']) is not None,
            f'{path}: version must be an ISO calendar date (YYYY-MM-DD)')
    try:
        date.fromisoformat(paper['version'])
    except ValueError as error:
        raise ValueError(f'{path}: invalid version date') from error
    for key in ('authors', 'topics'):
        value = paper.get(key)
        require(isinstance(value, list) and all(text(item) for item in value),
                f'{path}: {key} must be a list of nonempty strings')
    require(bool(paper['authors']), f'{path}: authors cannot be empty')
    doi = paper.get('doi')
    require(doi is None or text(doi), f'{path}: doi must be a nonempty string or null')
    if 'rights' in paper:
        require(isinstance(paper['rights'], dict), f'{path}: rights must be an object')
        for key, value in paper['rights'].items():
            require(value is None or isinstance(value, str),
                    f'{path}: rights.{key} must be a string or null')
    validation = paper.get('validation')
    require(isinstance(validation, dict) and text(validation.get('scope')) and
            type(validation.get('formal_verification')) is bool,
            f'{path}: validation needs scope (string) and formal_verification (boolean)')
    artifacts = paper.get('artifacts')
    require(isinstance(artifacts, list) and bool(artifacts),
            f'{path}: artifacts must be a nonempty list')
    listed = set()
    for artifact in artifacts:
        require(isinstance(artifact, dict), f'{path}: each artifact must be an object')
        relative = artifact.get('path')
        target = safe_file(path.parent, relative)
        require(relative not in listed, f'{path}: duplicate artifact: {relative}')
        listed.add(relative)
        require(text(artifact.get('role')), f'{path}: artifact role must be a nonempty string')
        require(isinstance(artifact.get('sha256'), str) and
                HASH.fullmatch(artifact['sha256']) is not None,
                f'{path}: artifact sha256 must contain 64 lowercase hexadecimal characters')
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        if refresh:
            artifact['sha256'] = digest
        else:
            require(digest == artifact['sha256'], f'{path}: checksum mismatch for {relative}')
    require('paper.pdf' in listed, f'{path}: paper.pdf must be listed as an artifact')
    require('source/README.md' in listed and any(
            name.startswith('source/') and name.endswith('.tex') for name in listed),
            f'{path}: list source/README.md and the author TeX source')
    safe_file(path.parent, 'README.md')
    # Build products are never adopted as artifacts. All other unlisted files
    # require an explicit inventory edit or removal, including during refresh.
    for target in path.parent.rglob('*'):
        relative = target.relative_to(path.parent)
        if 'build' in relative.parts or '__pycache__' in relative.parts:
            continue
        if target.is_symlink():
            raise ValueError(f'{path}: symlink in paper directory: {relative}')
        if target.is_file() and relative.as_posix() not in listed | WRAPPERS:
            if target.suffix in GENERATED_SUFFIXES:
                continue
            raise ValueError(f'{path}: unlisted file: {relative}; list it explicitly or remove it')
    return paper


def collect(root=ROOT, refresh_slug=None):
    if refresh_slug is not None:
        require(SLUG.fullmatch(refresh_slug) is not None, 'invalid refresh slug')
    paths = sorted((root / 'papers').glob('*/paper.json'))
    require(bool(paths), 'no paper metadata found')
    papers = [validate_paper(path, refresh=path.parent.name == refresh_slug) for path in paths]
    require(refresh_slug is None or any(p['slug'] == refresh_slug for p in papers),
            f'unknown paper slug: {refresh_slug}')
    require(len({p['slug'] for p in papers}) == len(papers), 'duplicate paper slugs')
    return papers


def refresh_hashes(root, slug):
    # Validate every paper before making the one explicitly requested write.
    papers = collect(root, refresh_slug=slug)
    paper = next(p for p in papers if p['slug'] == slug)
    path = root / 'papers' / slug / 'paper.json'
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                     prefix='.paper-', suffix='.tmp', delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(paper, stream, indent=2, ensure_ascii=False)
        stream.write('\n')
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)



def outputs(papers):
    lines = ['# Haibo Lu — Mathematical Research', '',
             'Papers, author sources, and reproducibility materials by Haibo Lu.', '',
             'Each paper has a stable directory with its author source and',
             'the scope of any supporting computation recorded in metadata.', '',
             '<!-- Generated by scripts/catalog.py from papers/*/paper.json. -->', '',
             '| Paper | Version | Materials |',
             '| --- | --- | --- |']
    for p in papers:
        base = 'papers/' + p['slug']
        title = p['title'].replace('|', r'\|')
        links = f'[PDF]({base}/paper.pdf) · [Source]({base}/source/) · [Details]({base}/README.md)'
        if p.get('doi'):
            links += f" · [DOI](https://doi.org/{p['doi']})"
        lines.append(f"| {title} | {p['version']} | {links} |")
    lines += ['', '## Maintenance', '',
              'See [CONTRIBUTING.md](CONTRIBUTING.md) to update a paper or add a new one.',
              'Run `python3 scripts/catalog.py` to regenerate the index and',
              '`python3 scripts/catalog.py --check` to validate it.', '',
              'Rights are stated per paper and artifact. This collection does not',
              'apply a new blanket license. Preserve author AI-use disclosures and',
              'distinguish analytic proofs, finite checks, and formal verification.', '']
    return {'README.md': '\n'.join(lines),
            'catalog.json': json.dumps({'schema_version': 1, 'papers': papers},
                                       indent=2, ensure_ascii=False) + '\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check', action='store_true', help='Fail if generated files are stale')
    mode.add_argument('--refresh-hashes', metavar='SLUG',
                      help='Update only listed artifact hashes for one paper; do not generate or publish')
    args = parser.parse_args()
    try:
        if args.refresh_hashes is not None:
            refresh_hashes(ROOT, args.refresh_hashes)
            print(f'Refreshed listed artifact hashes for {args.refresh_hashes}. '
                  'Run python3 scripts/catalog.py to regenerate indexes. Nothing was published.')
            return
        papers = collect()
        for name, content in outputs(papers).items():
            path = ROOT / name
            if args.check:
                require(path.exists() and path.read_text(encoding='utf-8') == content,
                        f'{name} is stale; run python3 scripts/catalog.py')
            else:
                path.write_text(content, encoding='utf-8')
        print(f'Validated {len(papers)} paper(s), all artifact hashes and generated indexes.')
    except (ValueError, OSError) as error:
        parser.exit(1, f'error: {error}\n')


if __name__ == '__main__':
    main()
