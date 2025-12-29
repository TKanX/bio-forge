/**
 * @file Home Page
 *
 * Landing page with hero, features, and links.
 */

import { HomeLayout } from "@/layouts";
import { Hero, Features, Links, Footer } from "@/features/home";

// ============================================================================
// Page
// ============================================================================

export default function HomePage() {
  return (
    <HomeLayout>
      <Hero />
      <Features />
      <Links />
      <Footer />
    </HomeLayout>
  );
}
