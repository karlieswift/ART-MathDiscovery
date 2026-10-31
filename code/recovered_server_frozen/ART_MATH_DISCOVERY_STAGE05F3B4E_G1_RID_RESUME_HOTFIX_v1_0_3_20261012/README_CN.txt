Stage05F3B4E G1-R3 existing-RID resume hotfix.

Purpose:
- R2 successfully submitted RT5 and obtained RID CR73H1ST014.
- The client then timed out after 10800 s while NCBI still reported the RID pending.
- R3 treats this as client-side transport waiting only, not a failed biological search.

Frozen scientific scope unchanged.
R3 changes only state handling:
1. reuse an existing valid RID instead of resubmitting it;
2. query SearchInfo once per minute;
3. allow up to 18 h additional waiting in one invocation;
4. on another client timeout, preserve the RID and run this same script again;
5. never interpret WAITING/timeout as no-hit.

NCBI RIDs are time-limited; do not intentionally leave a pending RID unattended for a full day.
