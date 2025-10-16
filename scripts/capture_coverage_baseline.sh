#!/bin/bash
# Capture coverage baseline for Phase 2 Testing & Security
# This script captures the current test coverage as a baseline to prevent regressions

set -e  # Exit on error

echo "📊 Capturing coverage baseline for Phase 2..."
echo ""

# Run full test suite with coverage
echo "Running test suite with coverage..."
pytest --cov=muc_one_up \
       --cov-report=json:coverage.baseline.json \
       --cov-report=term-missing \
       -v

echo ""
echo "✅ Coverage data captured to coverage.baseline.json"

# Check if jq is available for creating snapshot
if command -v jq &> /dev/null; then
    echo "Creating baseline snapshot..."

    # Create baseline snapshot with timestamp
    cat coverage.baseline.json | jq '{
      "timestamp": now | strftime("%Y-%m-%d %H:%M:%S"),
      "total_coverage": .totals.percent_covered,
      "files": [.files | to_entries[] | {
        "file": .key,
        "coverage": .value.summary.percent_covered,
        "missing_lines": .value.summary.missing_lines
      }]
    }' > coverage.baseline.snapshot.json

    echo "✅ Baseline snapshot created: coverage.baseline.snapshot.json"

    # Display baseline coverage
    total_coverage=$(jq '.totals.percent_covered' coverage.baseline.json)
    echo ""
    echo "📈 Baseline Coverage: ${total_coverage}%"
    echo ""
    echo "Next steps:"
    echo "  1. Commit coverage.baseline.json to git"
    echo "  2. Run 'python scripts/check_coverage_ratchet.py' after making changes"
    echo "  3. Coverage should not decrease below this baseline"
else
    echo "⚠️  jq not found - skipping snapshot creation"
    echo "   Install jq for enhanced baseline reporting: apt-get install jq"
    echo ""
    echo "✅ Baseline captured in coverage.baseline.json"
    echo ""
    echo "Next steps:"
    echo "  1. Commit coverage.baseline.json to git"
    echo "  2. Run 'python scripts/check_coverage_ratchet.py' after making changes"
fi

echo ""
echo "Done! 🎉"
