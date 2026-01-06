struct TraversalRule
    tags::Set{Pair{String, String}}
    traversable::Bool
end

TraversalRule(tags::AbstractVector{Pair{String, String}}, traversable::Bool) =
    TraversalRule(Set(tags), traversable)

struct RuleBasedTraversalPermissionSettings
    rules::Vector{TraversalRule}
    default::Bool
end

const DEFAULT_TRAVERSAL_SETTINGS = RuleBasedTraversalPermissionSettings([
    # things designated pedestrian or explicitly tagged as such - treat as walkable
    TraversalRule("foot" .=> ["yes", "designated"], true),

    # things designated not pedestrian
    TraversalRule(["foot", "access"] .=> "no", false),

    # things we assume are pedestrian
    TraversalRule([
            ("highway" .=> ["footway", "cycleway", "pedestrian", "track", "sidewalk", "service", "road", "steps", "path", "crossing", "residential"])...,
            ("sidewalk" .=> ["yes", "both", "left", "right"])...,
            "sidewalk:left" => "yes",
            "sidewalk:right" => "yes",
            "sidewalk:both" => "yes"
        ], true)
], false)

function is_traversable(settings::RuleBasedTraversalPermissionSettings, way)
    for rule in settings.rules
        for tag in pairs(way.tags)
            if tag ∈ rule.tags
                return rule.traversable
            end
        end
    end
    return settings.default
end