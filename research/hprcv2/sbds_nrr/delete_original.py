#!/usr/bin/env python3
"""Withdrawn unsafe deletion helper; use review_quarantine.py instead.

The original implementation tried to hash unreadable originals and renamed
quarantine directories while claiming freed storage. It did not implement
actual deletion or all trust gates. No deletion should use that code.
The reviewed replacement uses narrow regular-file unlink operations only after
independent replacement validation, COMPLETE audited/scored search checks,
supervisor success, and a pinned reviewed manifest. See review_quarantine.py.
"""
import sys

if __name__ == '__main__':
    sys.exit('REFUSED: use review_quarantine.py --manifest REVIEW.json, then '
             '--apply-sha256 REVIEW_DIGEST only after reviewing its live evidence.')
