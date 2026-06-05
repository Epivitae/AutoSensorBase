/**
 * AutoSensorBase — Cloudflare Worker for Visitor Tracking
 *
 * Tracks unique page visits using Cloudflare's request.cf geo-data
 * and stores counters + recent locations in KV.
 *
 * Deployment:
 *   1. npm create cloudflare@latest
 *   2. npx wrangler kv:namespace create VISITORS_KV
 *   3. npx wrangler kv:namespace create VISITORS_KV --preview  (for dev)
 *   4. Update wrangler.toml with the KV namespace ID
 *   5. npx wrangler deploy
 *
 * Local dev: npx wrangler dev  (request.cf is simulated in dev)
 */

export default {
  async fetch(request, env, ctx) {
    // ── CORS Preflight ──────────────────────────────────────────
    if (request.method === "OPTIONS") {
      return new Response(null, {
        status: 204,
        headers: corsHeaders(),
      });
    }

    // ── Only GET ────────────────────────────────────────────────
    if (request.method !== "GET") {
      return new Response("Method Not Allowed", {
        status: 405,
        headers: corsHeaders(),
      });
    }

    // ── Geo Data from Cloudflare ────────────────────────────────
    const cf = request.cf || {};
    const lat = cf.latitude;
    const lng = cf.longitude;

    // ── Bot Filtering (skip known crawlers) ─────────────────────
    const ua = (request.headers.get("user-agent") || "").toLowerCase();
    const isBot =
      ua.includes("bot") ||
      ua.includes("crawler") ||
      ua.includes("spider") ||
      ua.includes("scanner") ||
      ua.includes("headless") ||
      ua.includes("preview"); // Slack/Discord link previews

    // ── Increment Total Count ───────────────────────────────────
    let total = 0;
    try {
      const raw = await env.VISITORS_KV.get("total_count");
      total = parseInt(raw, 10) || 0;
      if (!isBot) {
        total += 1;
        await env.VISITORS_KV.put("total_count", total.toString());
      }
    } catch (e) {
      console.error("KV read/write error (total_count):", e);
    }

    // ── Append Location (if geo available & not a bot) ─────────
    let locations = [];
    try {
      const stored = await env.VISITORS_KV.get("recent_locations");
      locations = stored ? JSON.parse(stored) : [];

      if (!isBot && lat != null && lng != null) {
        // Round to ~1km precision (3 decimal places)
        locations.push({
          lat: Math.round(lat * 1000) / 1000,
          lng: Math.round(lng * 1000) / 1000,
        });

        // Keep only the most recent 100
        if (locations.length > 100) {
          locations = locations.slice(-100);
        }

        await env.VISITORS_KV.put("recent_locations", JSON.stringify(locations));
      }
    } catch (e) {
      console.error("KV read/write error (recent_locations):", e);
    }

    // ── Response ────────────────────────────────────────────────
    return new Response(
      JSON.stringify({
        total,
        locations,
      }),
      {
        status: 200,
        headers: {
          ...corsHeaders(),
          "Content-Type": "application/json",
          "Cache-Control": "no-cache, no-store, must-revalidate",
        },
      }
    );
  },
};

/**
 * CORS headers — allow your GitHub Pages domain.
 * Replace with your actual GitHub Pages URL after deployment.
 */
function corsHeaders() {
  return {
    "Access-Control-Allow-Origin": "*",
    "Access-Control-Allow-Methods": "GET, OPTIONS",
    "Access-Control-Allow-Headers": "Content-Type",
    "Access-Control-Max-Age": "86400",
  };
}
