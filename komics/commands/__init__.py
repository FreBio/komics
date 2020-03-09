from pkg_resources import get_distribution

try:
    __version__ = get_distribution('komics').version
except:
    __version__ = 'local'


from komics.commands import all, bam2fq, trimfq, assemble, circularize, polish, qualcheck