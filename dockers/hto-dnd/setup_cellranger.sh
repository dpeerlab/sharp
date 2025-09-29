PATH_CR=/Users/krauset/projects/sharp_dev/dockers/cellranger

# get download-url from
# https://www.10xgenomics.com/support/software/cell-ranger/downloads

# download to /tmp
cd /tmp
curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1749527681&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=Zbbk4lVAaXVmgJzjT4uAW2IRrkZOpIUG-gAbeoopZ-B3TM5QuZzh7nOlJUnyphfMB5ExSdNp9W7I3FWeBnlgbkiCYVlZMdnULl2CRbIy0BNWDNtnXJZJJ86G4XyvfwWLiEhadAoOwzKQjEHXhN3ozYwlQOzbUcJk4Nbf0xgrQW-BKr-LMDPpoT2GGyyXxPcS2qVuQVpe1cO2q3rAphLLWwr48AtXgvdSDRl-a9D7MkV-dnTJfAJdVB-1Ql4x6hwW26ld6febF8tPDdN9WqMJHUHngNRZJA0U5tYw1lir4ka4VOb74mMg8KiKXRRnKENVa9FZwUMjTz88YW94zVltQg__"

# extract to $PATH_CR
tar -xzvf cellranger-9.0.1.tar.gz -C $PATH_CR