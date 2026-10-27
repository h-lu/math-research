"""Exercise safe catalog maintenance against disposable paper directories."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'catalog.py'


class CatalogTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'scripts').mkdir()
        shutil.copyfile(SCRIPT, self.root / 'scripts/catalog.py')
        self.paper('example-paper')
        self.run_cli()

    def paper(self, slug):
        base = self.root / 'papers' / slug
        (base / 'source').mkdir(parents=True)
        files = {'paper.pdf': b'%PDF-1.7\nfixture', 'source/main.tex': b'Author source',
                 'source/README.md': b'Build instructions', 'README.md': b'Paper guide'}
        for name, content in files.items():
            (base / name).write_bytes(content)
        data = {'schema_version': 1, 'slug': slug, 'title': 'Example title',
                'authors': ['Example Author'], 'version': '2026-10-27', 'abstract': 'Example.',
                'topics': [], 'doi': None,
                'validation': {'scope': 'Fixture only', 'formal_verification': False},
                'artifacts': [{'path': name, 'role': 'author-source',
                               'sha256': hashlib.sha256(content).hexdigest()}
                              for name, content in files.items() if name != 'README.md']}
        (base / 'paper.json').write_text(json.dumps(data))
        return base

    def run_cli(self, *args, success=True):
        result = subprocess.run([sys.executable, 'scripts/catalog.py', *args], cwd=self.root,
                                text=True, capture_output=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def change_metadata(self, change):
        path = self.root / 'papers/example-paper/paper.json'
        data = json.loads(path.read_text())
        change(data)
        path.write_text(json.dumps(data))

    def test_missing_file(self):
        (self.root / 'papers/example-paper/paper.pdf').unlink()
        self.assertIn('missing required file', self.run_cli('--check', success=False).stderr)

    def test_wrong_hash_is_not_silently_updated(self):
        path = self.root / 'papers/example-paper/paper.json'
        before = path.read_bytes()
        (path.parent / 'source/main.tex').write_text('Changed source')
        self.assertIn('checksum mismatch', self.run_cli(success=False).stderr)
        self.assertEqual(path.read_bytes(), before)

    def test_parent_escape(self):
        self.change_metadata(lambda d: d['artifacts'][0].update(path='../outside.pdf'))
        self.assertIn('unsafe artifact path', self.run_cli('--refresh-hashes', 'example-paper', success=False).stderr)

    def test_absolute_path(self):
        self.change_metadata(lambda d: d['artifacts'][0].update(path=str(self.root / 'outside.pdf')))
        self.assertIn('unsafe artifact path', self.run_cli('--check', success=False).stderr)

    def test_symlink_ancestor(self):
        base = self.root / 'papers/example-paper'
        shutil.move(base / 'source', self.root / 'external-source')
        (base / 'source').symlink_to(self.root / 'external-source', target_is_directory=True)
        self.assertIn('artifact path', self.run_cli('--check', success=False).stderr)

    def test_stale_index(self):
        (self.root / 'README.md').write_text('stale')
        self.assertIn('stale', self.run_cli('--check', success=False).stderr)

    def test_refresh_only_selected_paper_and_no_index_generation(self):
        second = self.paper('another-paper')
        second_before = (second / 'paper.json').read_bytes()
        index_before = (self.root / 'catalog.json').read_bytes()
        path = self.root / 'papers/example-paper/source/main.tex'
        path.write_text('Revised author source')
        self.run_cli('--refresh-hashes', 'example-paper')
        data = json.loads((path.parent.parent / 'paper.json').read_text())
        artifact = next(a for a in data['artifacts'] if a['path'] == 'source/main.tex')
        self.assertEqual(artifact['sha256'], hashlib.sha256(path.read_bytes()).hexdigest())
        self.assertEqual((second / 'paper.json').read_bytes(), second_before)
        self.assertEqual((self.root / 'catalog.json').read_bytes(), index_before)
        self.run_cli('--check', success=False)
        self.run_cli()
        self.run_cli('--check')

    def test_refresh_refuses_unlisted_file_without_writing(self):
        base = self.root / 'papers/example-paper'
        before = (base / 'paper.json').read_bytes()
        (base / 'source/unlisted.tex').write_text('Not approved for inventory')
        self.assertIn('unlisted file', self.run_cli('--refresh-hashes', 'example-paper', success=False).stderr)
        self.assertEqual((base / 'paper.json').read_bytes(), before)

    def test_required_guides(self):
        (self.root / 'papers/example-paper/README.md').unlink()
        self.run_cli('--check', success=False)

    def test_basic_types_and_calendar_dates(self):
        path = self.root / 'papers/example-paper/paper.json'
        original = path.read_bytes()
        for key, value in [('authors', 'Author'), ('topics', [1]), ('version', '2026-02-30'),
                           ('slug', '../example-paper'), ('schema_version', True),
                           ('validation', {'scope': 'Checks', 'formal_verification': 'no'})]:
            with self.subTest(key=key):
                path.write_bytes(original)
                self.change_metadata(lambda d: d.update({key: value}))
                self.run_cli('--check', success=False)

    def test_unknown_refresh_slug(self):
        self.run_cli('--refresh-hashes', 'missing-paper', success=False)

    def test_modes_are_mutually_exclusive(self):
        self.run_cli('--check', '--refresh-hashes', 'example-paper', success=False)


if __name__ == '__main__':
    unittest.main()
