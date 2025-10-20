#!/bin/bash
# Bulk dismiss Trivy code scanning alerts
# These are base image vulnerabilities that cannot be fixed

REPO="berntpopp/MucOneUp"
REASON="won't fix"
COMMENT="Base image vulnerability with no available fix. Filtered by ignore-unfixed flag in workflow."

echo "Fetching all open code scanning alerts (this may take a moment)..."

# Fetch ALL alerts with proper pagination using jq to process per-page
ALERT_NUMBERS=$(gh api \
    --paginate \
    -X GET \
    repos/$REPO/code-scanning/alerts \
    -f state=open \
    -f tool_name=Trivy \
    -f per_page=100 \
    --jq '.[] | .number')

TOTAL=$(echo "$ALERT_NUMBERS" | wc -l | tr -d ' ')
echo "Found $TOTAL open Trivy alerts to dismiss"

if [ "$TOTAL" -eq 0 ]; then
    echo "No alerts to dismiss"
    exit 0
fi

echo "Starting bulk dismissal..."
echo ""

COUNT=0
FAILED=0

for ALERT_NUM in $ALERT_NUMBERS; do
    COUNT=$((COUNT + 1))

    # Show progress every 50 alerts
    if [ $((COUNT % 50)) -eq 0 ] || [ $COUNT -eq 1 ]; then
        echo "Progress: $COUNT/$TOTAL alerts processed..."
    fi

    gh api \
        --method PATCH \
        -H "Accept: application/vnd.github+json" \
        repos/$REPO/code-scanning/alerts/$ALERT_NUM \
        -f state='dismissed' \
        -f dismissed_reason="$REASON" \
        -f dismissed_comment="$COMMENT" \
        --silent > /dev/null 2>&1

    if [ $? -ne 0 ]; then
        FAILED=$((FAILED + 1))
    fi

    # Rate limiting: brief pause every 100 requests
    if [ $((COUNT % 100)) -eq 0 ]; then
        echo "  (pausing 5s for rate limiting...)"
        sleep 5
    fi
done

echo ""
echo "================================"
echo "✓ Dismissed: $((COUNT - FAILED)) alerts"
if [ $FAILED -gt 0 ]; then
    echo "✗ Failed: $FAILED alerts"
fi
echo "================================"
