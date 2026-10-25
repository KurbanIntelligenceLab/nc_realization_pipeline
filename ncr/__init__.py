"""ncr: data pipeline for "Computational predictions of stable inorganic materials and
their later experimental synthesis". Every step reads files you provide and writes
files under the work directory; nothing is downloaded except, if you enable it, the
literature search (OpenAlex) and the language-model calls (Anthropic API)."""
__version__ = "2.2.2"
